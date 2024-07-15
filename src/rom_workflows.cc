// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "rom_workflows.h"
#include "Electrostatic.h"
#include <memory>
#include <string>
#include <stdexcept>

template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    auto size = static_cast<size_t>( size_s );
    std::unique_ptr<char[]> buf( new char[ size ] );
    std::snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

template <class OrbitalsType>
void readRestartFiles(MGmolInterface *mgmol_)
{
    Control& ct              = *(Control::instance());
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();

    ROMPrivateOptions rom_options = ct.getROMOptions();
    /* type of variable we intend to run POD */
    ROMVariable rom_var = rom_options.variable;

    /* number of restart files, start/end indices */
    assert(rom_options.restart_file_minidx >= 0);
    assert(rom_options.restart_file_maxidx >= 0);
    const int minidx = rom_options.restart_file_minidx;
    const int maxidx = rom_options.restart_file_maxidx;
    const int num_restart = maxidx - minidx + 1;

    MGmol<OrbitalsType> *mgmol = static_cast<MGmol<OrbitalsType> *>(mgmol_);
    OrbitalsType *orbitals = mgmol->getOrbitals();
    Potentials& pot = mgmol->getHamiltonian()->potential();
    std::string filename;

    /* Determine basis prefix, dimension, and sample size */
    std::string basis_prefix = rom_options.basis_file;
    int dim;
    int totalSamples = num_restart;
    const int chrom_num = orbitals->chromatic_number();
    switch (rom_var)
    {
    case ROMVariable::ORBITALS:
        basis_prefix += "_orbitals";
        dim = orbitals->getLocNumpt();
        /* if orbitals, each sample have chromatic number of wave functions */
        totalSamples *= orbitals->chromatic_number();
        break;

    case ROMVariable::POTENTIAL:
        basis_prefix += "_potential";
        dim = pot.size();
        break;
    
    default:
        ct.global_exit(-1);
        break;
    }

    /* Initialize libROM classes */
    CAROM::Options svd_options(dim, totalSamples, 1);
    CAROM::BasisGenerator basis_generator(svd_options, false, basis_prefix);
    
    /* Collect the restart files */
    for (int k = minidx; k <= maxidx; k++)
    {
        filename = string_format(rom_options.restart_file_fmt, k);
        mgmol->loadRestartFile(filename);
        assert(dim == orbitals->getLocNumpt());
        assert(chrom_num == orbitals->chromatic_number());

        switch (rom_var)
        {
        case ROMVariable::ORBITALS:
            for (int i = 0; i < chrom_num; ++i)
                basis_generator.takeSample(orbitals->getPsi(i));
            break;

        case ROMVariable::POTENTIAL:
            basis_prefix += "_potential";
            /* we save hartree potential */
            basis_generator.takeSample(pot.vh_rho());
            break;
        }
    }
    basis_generator.writeSnapshot();
    basis_generator.endSamples();
}

template <class OrbitalsType>
void buildROMPoissonOperator(MGmolInterface *mgmol_)
{
    Control& ct              = *(Control::instance());
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();

    ROMPrivateOptions rom_options = ct.getROMOptions();
    /* type of variable we intend to run POD */
    ROMVariable rom_var = rom_options.variable;
    if (rom_var != ROMVariable::POTENTIAL)
    {
        std::cerr << "buildROMPoissonOperator error: ROM variable must be POTENTIAL to run this stage!\n" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    /* Load Hartree potential basis matrix */
    std::string basis_file = rom_options.basis_file;
    const int num_pot_basis = rom_options.num_potbasis;
    CAROM::BasisReader basis_reader(basis_file);
    CAROM::Matrix *pot_basis = basis_reader.getSpatialBasis(num_pot_basis);

    /* Load PoissonSolver pointer */
    MGmol<OrbitalsType> *mgmol = static_cast<MGmol<OrbitalsType> *>(mgmol_);
    Poisson *poisson = mgmol->electrostat_->getPoissonSolver();

    /* Initialize ROM matrix (undistributed) */
    CAROM::Matrix pot_rom(num_pot_basis, num_pot_basis, false);

    pb::GridFunc<POTDTYPE> gf_col(poisson->vh());
    pb::GridFunc<POTDTYPE> gf_opcol(gf_col);
    CAROM::Vector *col = nullptr, *op_col = nullptr, *rom_col = nullptr;
    for (int c = 0; c < num_pot_basis; c++)
    {
        /* copy c-th column librom vector to GridFunc gf_col */
        col = pot_basis->getColumn(c);
        gf_col.assign(col->getData());

        /* apply Laplace operator */
        poisson->applyOperator(gf_col, gf_opcol);

        /* get librom view-vector of gf_opcol */
        op_col = new CAROM::Vector(gf_opcol.uu(), col->dim(), true, false);

        /* Compute basis projection of the column */
        /* Resulting vector is undistributed */
        rom_col = pot_basis->transposeMult(*op_col);

        /* libROM matrix is row-major, so data copy is necessary */
        for (int r = 0; r < num_pot_basis; r++)
            pot_rom(r, c) = (*rom_col)(r);

        delete col;
        delete op_col;
        delete rom_col;
    }   // for (int c = 0; c < num_pot_basis; c++)
    
    /* Save ROM operator */
    // write the file from PE0 only
    if (MPIdata::onpe0)
    {
        std::string rom_oper = "pot_rom_oper.h5";
        CAROM::HDFDatabase h5_helper;
        h5_helper.create(rom_oper);
        h5_helper.putInteger("number_of_potential_basis", num_pot_basis);
        h5_helper.putDoubleArray("potential_rom_operator", pot_rom.getData(),
                                num_pot_basis * num_pot_basis, false);

        /* save the inverse as well */
        pot_rom.inverse();
        h5_helper.putDoubleArray("potential_rom_inverse", pot_rom.getData(),
                                num_pot_basis * num_pot_basis, false);

        h5_helper.close();
    }
}

template <class OrbitalsType>
void testROMPoissonOperator(MGmolInterface *mgmol_)
{
    Control& ct              = *(Control::instance());
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();

    ROMPrivateOptions rom_options = ct.getROMOptions();

    /* Load MGmol pointer and Potentials */
    MGmol<OrbitalsType> *mgmol = static_cast<MGmol<OrbitalsType> *>(mgmol_);
    Poisson *poisson = mgmol->electrostat_->getPoissonSolver();
    Potentials& pot = mgmol->getHamiltonian()->potential();
    const int dim = pot.size();
    printf("pot size: %d\n", dim);

    /* GridFunc initialization inputs */
    const pb::Grid &grid(poisson->vh().grid());
    short bc[3];
    for (int d = 0; d < 3; d++)
        bc[d] = poisson->vh().bc(d);

    /* fictitious snapshot numbers */
    const int nsnapshot = 3;

    /* Set compensating charges to zero now */
    pb::GridFunc<POTDTYPE> rhoc(grid, bc[0], bc[1], bc[2]);
    rhoc = 0.0;

    /* Generate fictitious right-hand sides and snapshots */
    std::vector<std::vector<POTDTYPE>> rhs(nsnapshot), fom_sol(nsnapshot);
    for (int s = 0; s < nsnapshot; s++)
    {
        rhs[s].resize(dim);
        for (int d = 0; d < dim; d++)
            rhs[s][d] = ran0();
        
        /* average out for periodic bc */
        pb::GridFunc<POTDTYPE> rhs_gf(grid, bc[0], bc[1], bc[2]);
        rhs_gf.assign(rhs[s].data());
        double avg = rhs_gf.get_average();
        rhs_gf -= avg;

        /* copy back to rhs */
        rhs_gf.init_vect(rhs[s].data(), 'd');

        poisson->solve(rhs_gf, rhoc);

        fom_sol[s].resize(dim);
        poisson->vh().init_vect(fom_sol[s].data(), 'd');

        /* check if the solution is correct */
        pb::GridFunc<POTDTYPE> res(grid, bc[0], bc[1], bc[2]);
        pb::GridFunc<POTDTYPE> sol_gf(grid, bc[0], bc[1], bc[2]);
        sol_gf.assign(fom_sol[s].data());
        /* apply Laplace operator */
        poisson->applyOperator(sol_gf, res);
        /* FD operator scales rhs by 4pi */
        res.axpy(- 4. * M_PI, rhs_gf);
        printf("FOM res norm: %.3e\n", res.norm2());
    }

    // /* Initialize libROM classes */
    // std::string basis_prefix = "test_poisson";
    // CAROM::Options svd_options(dim, nsnapshot, 1);
    // CAROM::BasisGenerator basis_generator(svd_options, false, basis_prefix);

    // /* Collect snapshots and train POD basis */
    // for (int s = 0; s < nsnapshot; s++)
    //     basis_generator.takeSample(fom_sol[s]->uu());
    // basis_generator.endSamples();

    // /* Load POD basis. We use the maximum number of basis vectors. */
    // const CAROM::Matrix *pot_basis = basis_generator.getSpatialBasis();

    // /* Check if full projection preserves FOM solution */
    // for (int c = 0; c < nsnapshot; c++)
    // {
    //     CAROM::Vector *fom_sol_vec = nullptr;
    //     /* get librom view-vector of fom_sol */
    //     fom_sol_vec = new CAROM::Vector(fom_sol[c]->uu(), pot_basis->numRows(), true, false);

    //     CAROM::Vector *rom_proj = pot_basis->transposeMult(*fom_sol_vec);
    //     CAROM::Vector *reconstruct = pot_basis->mult(*rom_proj);

    //     /* error on libROM side */
    //     CAROM::Vector *librom_error = reconstruct->minus(fom_sol_vec);
    //     printf("librom reconstruction error: %.3e\n", librom_error->norm());

    //     /* error on mgmol side */
    //     pb::GridFunc<POTDTYPE> gf_col(grid, bc[0], bc[1], bc[2]);
    //     printf("reconstruct size: %d\n", reconstruct->dim());
    //     printf("gf_col sizeg: %d\n", gf_col.sizeg());
    //     printf("gf_col size: %d\n", gf_col.size());
    //     gf_col.setValues(reconstruct->dim(), reconstruct->getData());
    //     gf_col -= *fom_sol[c];
    //     printf("mgmol reconstruction error: %.3e\n", gf_col.norm2());

    //     delete fom_sol_vec;
    //     delete rom_proj;
    //     delete reconstruct;
    //     delete librom_error;
    // }

    // /* Check FOM axpy is equivalent to ROM axpy */
    // for (int s = 0; s < nsnapshot; s++)
    // {
    //     pb::GridFunc<POTDTYPE> res(grid, bc[0], bc[1], bc[2]);
    //     /* apply Laplace operator */
    //     poisson->applyOperator(*fom_sol[s], res);

    //     /* get librom view-vector of fom_res */
    //     CAROM::Vector *fom_res = new CAROM::Vector(res.uu(), pot_basis->numRows(), true, false);
    //     CAROM::Vector *rom_res = pot_basis->transposeMult(*fom_res);

    //     /* get librom view-vector of fom_rhs */
    //     pb::GridFunc<POTDTYPE> mgmol_rhs(*rhs[s]);
    //     CAROM::Vector *fom_rhs = new CAROM::Vector(mgmol_rhs.uu(), pot_basis->numRows(), true, false);
    //     CAROM::Vector *rom_rhs = pot_basis->transposeMult(*fom_rhs);

    //     /* FD operator scales rhs by 4pi */
    //     res.axpy(- 4. * M_PI, mgmol_rhs);
    //     printf("FOM res norm: %.3e\n", res.norm2());
    //     delete fom_res;
    //     fom_res = new CAROM::Vector(res.uu(), pot_basis->numRows(), true, false);
    //     CAROM::Vector *res_proj = pot_basis->transposeMult(*fom_res);
    //     printf("FOM res projection norm: %.3e\n", res_proj->norm());

    //     *rom_rhs *= 4. * M_PI;
    //     *rom_res -= *rom_rhs;

    //     printf("ROM res norm: %.3e\n", rom_res->norm());

    //     delete fom_res;
    //     delete rom_res;
    //     delete fom_rhs;
    //     delete rom_rhs;
    // }

    // /* Initialize Projection ROM matrix (undistributed) */
    // CAROM::Matrix pot_rom(nsnapshot, nsnapshot, false);

    // /* Build Projection of Poisson operator */
    // pb::GridFunc<POTDTYPE> gf_col(poisson->vh());
    // pb::GridFunc<POTDTYPE> gf_opcol(gf_col);
    // CAROM::Vector *col = nullptr, *op_col = nullptr, *rom_col = nullptr;
    // for (int c = 0; c < nsnapshot; c++)
    // {
    //     /* copy c-th column librom vector to GridFunc gf_col */
    //     col = pot_basis->getColumn(c);
    //     // gf_col.assign(col->getData());
    //     gf_col.setValues(col->dim(), col->getData());

    //     /* apply Laplace operator */
    //     poisson->applyOperator(gf_col, gf_opcol);

    //     /* get librom view-vector of gf_opcol */
    //     op_col = new CAROM::Vector(gf_opcol.uu(), col->dim(), true, false);

    //     /* Compute basis projection of the column */
    //     /* Resulting vector is undistributed */
    //     rom_col = pot_basis->transposeMult(*op_col);

    //     /* libROM matrix is row-major, so data copy is necessary */
    //     for (int r = 0; r < nsnapshot; r++)
    //         pot_rom(r, c) = (*rom_col)(r);

    //     delete col;
    //     delete op_col;
    //     delete rom_col;
    // }   // for (int c = 0; c < num_pot_basis; c++)

    // /* Inverse of the projection ROM matrix */
    // CAROM::Matrix pot_rom_inv(pot_rom);
    // pot_rom_inv.inverse();

    // /* Check the inverse */
    // CAROM::Matrix *identity = pot_rom_inv.mult(pot_rom);
    // printf("pot_rom_inv * pot_rom = identity\n");
    // for (int i = 0; i < nsnapshot; i++)
    // {
    //     for (int j = 0; j < nsnapshot; j++)
    //         printf("%.3e\t", identity->item(i, j));
    //     printf("\n");
    // }
    // delete identity;

    // /* Test with sample RHS. ROM must be able to 100% reproduce the FOM solution. */
    // std::vector<CAROM::Vector *> rom_sol(0), rom_rhs(0);
    // std::vector<pb::GridFunc<POTDTYPE> *> test_sol(0);
    // for (int s = 0; s < nsnapshot; s++)
    // {
    //     /* get librom view-vector of rhs[s] */
    //     op_col = new CAROM::Vector(rhs[s]->uu(), dim, true, false);

    //     /* project onto POD basis */
    //     rom_rhs.push_back(pot_basis->transposeMult(*op_col));
    //     delete op_col;

    //     /* FOM FD operator scales rhs by 4pi */
    //     *rom_rhs.back() *= 4. * M_PI;

    //     /* solve ROM */
    //     rom_sol.push_back(pot_rom_inv.mult(*rom_rhs.back()));

    //     /* check ROM solution */
    //     CAROM::Vector &res(*pot_rom.mult(*rom_sol.back()));
    //     res -= *rom_rhs.back();
    //     printf("rom res norm: %.3e\n", res.norm());

    //     /* initialize lift-up FOM solution */
    //     test_sol.push_back(new pb::GridFunc<POTDTYPE>(poisson->vh()));
    //     test_sol.back()->init_rand();

    //     /* lift back to FOM to check the accuracy */
    //     /* get librom view-vector of test_sol[s] */
    //     op_col = new CAROM::Vector(test_sol.back()->uu(), dim, true, false);
    //     pot_basis->mult(*rom_sol.back(), *op_col);

    //     delete op_col;
    // }

    // /* Compute relative errors */
    // for (int s = 0; s < nsnapshot; s++)
    // {
    //     *test_sol[s] -= *fom_sol[s];
    //     double rel_error = test_sol[s]->norm2() / fom_sol[s]->norm2();
    //     printf("%d-th sample relative error: %.3e\n", s, rel_error);
    // }

    /* clean up pointers */
    for (int s = 0; s < nsnapshot; s++)
    {
        // delete rom_sol[s];
        // delete rom_rhs[s];
        // delete test_sol[s];
    }
}

template void readRestartFiles<LocGridOrbitals>(MGmolInterface *mgmol_);
template void readRestartFiles<ExtendedGridOrbitals>(MGmolInterface *mgmol_);

template void buildROMPoissonOperator<LocGridOrbitals>(MGmolInterface *mgmol_);
template void buildROMPoissonOperator<ExtendedGridOrbitals>(MGmolInterface *mgmol_);

template void testROMPoissonOperator<LocGridOrbitals>(MGmolInterface *mgmol_);
template void testROMPoissonOperator<ExtendedGridOrbitals>(MGmolInterface *mgmol_);