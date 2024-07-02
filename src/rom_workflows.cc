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
            /* we save total potential for now */
            basis_generator.takeSample(pot.vtot());
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

template void readRestartFiles<LocGridOrbitals>(MGmolInterface *mgmol_);
template void readRestartFiles<ExtendedGridOrbitals>(MGmolInterface *mgmol_);

template void buildROMPoissonOperator<LocGridOrbitals>(MGmolInterface *mgmol_);
template void buildROMPoissonOperator<ExtendedGridOrbitals>(MGmolInterface *mgmol_);