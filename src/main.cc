// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

//
//                  main.cc
//
//    Description:
//        Real grid, finite difference, molecular dynamics program
//        for with nonorthogonal localized orbitals.
//
//        Uses Mehrstellen operators, multigrid accelerations, and
//        non-local pseudopotentials.
//
//     Includes LDA and PBE exchange and correlation functionals.
//
// Units:
//   Potentials, eigenvalues and operators in Rydberg
//   Energies in Hartree
//
#include <cassert>
#include <iostream>
#include <iterator>
#include <vector>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_CNR
#include <mkl.h>
#endif

#include <mpi.h>

#include "Control.h"
#include "DistMatrix.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "MatricesBlacsContext.h"
#include "Mesh.h"
#include "PackedCommunicationBuffer.h"
#include "ReplicatedWorkSpace.h"
#include "SparseDistMatrix.h"
#include "magma_singleton.h"
#include "tools.h"

#include <fenv.h>
#include <sys/cdefs.h>
#include <time.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

//#include "MemTrack.h"

/*
void trapfpe () {
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
*/

// A helper function
template <class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(cout, " "));
    return os;
}

int main(int argc, char** argv)
{
    // change handling of memory allocation errors
    set_new_handler(noMoreMemory);

    cout.sync_with_stdio();

    int mpirc = MPI_Init(&argc, &argv);
    if (mpirc != MPI_SUCCESS)
    {
        cerr << "MPI Initialization failed!!!" << endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    assert(mype > -1);
    onpe0 = (mype == 0);

#ifdef HAVE_MAGMA
    magma_int_t magmalog;

    magmalog = magma_init();
    if (magmalog == MAGMA_SUCCESS)
    {
        std::cout << "MAGMA Initialization: success" << std::endl;
    }
    else
    {
        if (magmalog == MAGMA_ERR_UNKNOWN)
            std::cout << "MAGMA Initialization: unknown error" << std::endl;
        if (magmalog == MAGMA_ERR_HOST_ALLOC)
            std::cout << "MAGMA Initialization: fails host alloc" << std::endl;
        return 1;
    }
#endif

    string input_file("");
    string lrs_filename;
    string constraints_filename("");
    bool tcheck = false;

    float total_spin = 0.;
    bool with_spin   = false;

    po::variables_map vm;

    // use configure file if it can be found
    // std::string config_filename("mgmol.cfg");

    // read options from PE0 only
    if (onpe0) try
        {
            string config_file;

            // Declare a group of options that will be
            // allowed only on command line
            po::options_description generic("Generic options");
            generic.add_options()("version,v", "print version string")("help,h",
                "produce help message")("check", "check input")("config,c",
                po::value<string>(&config_file)->default_value("mgmol.cfg"),
                "name of a file of a configuration.")("atomicCoordinates,i",
                po::value<vector<string>>(),
                "coordinates filename")("LRsFilename,l",
                po::value<string>(&lrs_filename), "LRs filename");

            // Declare a group of options (with default when appropriate) that
            // will be allowed in config file
            po::options_description config("Configuration");
            config.add_options()(
                "spin", po::value<float>(), "system total spin")("verbosity",
                po::value<short>()->default_value(1), "verbosity level")(
                "constraintsFilename", po::value<string>(&constraints_filename),
                "Name of file with list of constraints")("xcFunctional",
                po::value<string>()->required(), "XC functional: LDA or PBE")(
                "FDtype", po::value<string>()->default_value("Mehrstellen"),
                "Finite Difference scheme")("charge",
                po::value<short>()->default_value(0), "system total charge")(
                "Domain.lx", po::value<float>()->required(),
                "domain dimension in x direction")("Domain.ly",
                po::value<float>()->required(),
                "domain dimension in y direction")("Domain.lz",
                po::value<float>()->required(),
                "domain dimension in z direction")("Domain.ox",
                po::value<float>()->required(), "domain origin in x direction")(
                "Domain.oy", po::value<float>()->required(),
                "domain origin in y direction")("Domain.oz",
                po::value<float>()->required(), "domain origin in z direction")(
                "Mesh.nx", po::value<short>()->required(),
                "mesh dimension in x direction")("Mesh.ny",
                po::value<short>()->required(),
                "mesh dimension in y direction")("Mesh.nz",
                po::value<short>()->required(),
                "mesh dimension in z direction")("Potentials.pseudopotential",
                po::value<vector<string>>()->multitoken(),
                "pseudopotentials list")("Potentials.external",
                po::value<vector<string>>()->multitoken(),
                "external potentials list")("Potentials.binExternal",
                po::value<bool>()->default_value(true),
                "binary external potential")("Restart.input_filename",
                po::value<string>()->default_value(""),
                "Read restart filename/directory")("Restart.input_level",
                po::value<short>()->default_value(0),
                "Read restart level")("Restart.input_type",
                po::value<string>()->default_value("distributed"),
                "Read restart type: distributed or single_file")(
                "Restart.output_filename",
                po::value<string>()->default_value("auto"),
                "Dump restart filename/directory")("Restart.output_level",
                po::value<short>()->default_value(3),
                "Write restart level")("Restart.output_type",
                po::value<string>()->default_value("distributed"),
                "Write restart type: distributed or single_file")(
                "Restart.interval", po::value<short>()->default_value(1000),
                "Restart frequency")("Restart.rescale_v",
                po::value<double>()->default_value(1.),
                "rescaling factor velocity of all atoms")("Poisson.bcx",
                po::value<string>()->default_value("periodic"),
                "boundary condition x")("Poisson.bcy",
                po::value<string>()->default_value("periodic"),
                "boundary condition y")("Poisson.bcz",
                po::value<string>()->default_value("periodic"),
                "boundary condition z")("Poisson.diel",
                po::value<string>()->default_value("off"),
                "continuum solvent: on/off")("Run.type",
                po::value<string>()->default_value("QUENCH"), "Run type")(
                "Quench.solver", po::value<string>()->default_value("ABPG"),
                "Iterative solver for quench")("Quench.max_steps",
                po::value<short>()->default_value(200),
                "Max. steps in loose quench")("Quench.max_steps_tight",
                po::value<short>()->default_value(1000),
                "Max. steps in tight quench")("Quench.atol",
                po::value<float>()->default_value(1.e-12),
                "Abs. tol. in quench convergence")("Quench.rtol",
                po::value<float>()->default_value(-1.),
                "Rel. tol. in quench convergence")("Quench.conv_criterion",
                po::value<string>()->default_value("deltaE"),
                "Convergence criterion")("Quench.MLWC", po::value<bool>(),
                "Compute MLWC in quench")("Quench.MLWF",
                po::value<bool>()->default_value(false),
                "Compute MLWF (apply rotation) in quench")(
                "Quench.num_lin_iterations",
                po::value<short>()->default_value(0),
                "Number of iterations without potential update in quench")(
                "Quench.preconditioner_num_levels",
                po::value<short>()->default_value(2),
                "Number of levels for MG preconditioner")(
                "Quench.spread_penalty_damping",
                po::value<float>()->default_value(0.),
                "Spread penalty damping factor")("Quench.spread_penalty_target",
                po::value<float>()->default_value(2.),
                "Spread penalty target")("SpreadPenalty.type",
                po::value<string>()->default_value("individual"),
                "Spread penalty type (individual,volume)")(
                "SpreadPenalty.damping", po::value<float>()->default_value(1.),
                "Spread penalty damping factor")("SpreadPenalty.target",
                po::value<float>()->default_value(-1.),
                "Spread penalty target")("SpreadPenalty.alpha",
                po::value<float>()->default_value(0.),
                "Spread penalty factor")("MD.num_steps",
                po::value<short>()->default_value(1), "number of MD steps")(
                "MD.dt", po::value<float>(), "time step for MD (a.u.)")(
                "MD.print_interval", po::value<short>()->default_value(1),
                "print intervale for MD data")("MD.print_directory",
                po::value<string>()->default_value("MD"),
                "print directory for MD data")("MD.thermostat",
                po::value<string>()->default_value("OFF"),
                "MD thermostat: ON or OFF")("MD.remove_mass_center_motion",
                po::value<bool>()->default_value(true),
                "Remove mass center motion")("MD.type",
                po::value<string>()->default_value("BOMD"), "MD type: BOMD")(
                "GeomOpt.type", po::value<string>()->default_value("LBFGS"),
                "Geometry optimization algorithm")("GeomOpt.tol",
                po::value<float>()->default_value(4.e-4),
                "Tolerance on forces for Geometry optimization")(
                "GeomOpt.max_steps", po::value<short>()->default_value(1),
                "max. number of Geometry optimization steps")("GeomOpt.dt",
                po::value<float>(), "Delta t for trial pseudo-time steps")(
                "atomicCoordinates", po::value<vector<string>>(),
                "coordinates filename")("Thermostat.type",
                po::value<string>()->default_value("Langevin"),
                "Thermostat type")("Thermostat.temperature",
                po::value<float>()->default_value(-1.),
                "Thermostat temperature")("Thermostat.relax_time",
                po::value<float>()->default_value(-1.),
                "Thermostat relaxation time")("Thermostat.width",
                po::value<float>()->default_value(-1.),
                "Thermostat width (for SCALING)")("Orbitals.nempty",
                po::value<short>()->default_value(0),
                "Number of empty orbitals")("Orbitals.initial_type",
                po::value<string>()->default_value("random"),
                "initial orbitals type")("Orbitals.initial_width",
                po::value<float>()->default_value(10000.),
                "initial orbitals radius")("Orbitals.temperature",
                po::value<float>()->default_value(0.),
                "electronic temperature [K]")("ProjectedMatrices.solver",
                po::value<string>()->default_value("exact"),
                "solver for projected matrices")("ProjectedMatrices.printMM",
                po::value<bool>()->default_value(false),
                "print projected matrices in MM format")(
                "LocalizationRegions.radius",
                po::value<float>()->default_value(1000.),
                "Localization regions radius")("LocalizationRegions.adaptive",
                po::value<bool>()->default_value(true),
                "Localization regions adaptivity")(
                "LocalizationRegions.move_tol",
                po::value<float>()->default_value(1000.),
                "Localization regions move tolerance")(
                "Parallel.atomic_info_radius",
                po::value<float>()->default_value(8.),
                "Max. distance for atomic data communication")(
                "LoadBalancing.alpha", po::value<float>()->default_value(0.0),
                "Parameter for computing bias for load balancing algo")(
                "LoadBalancing.damping_tol",
                po::value<float>()->default_value(0.9),
                "Damping parameter for computing bias for load balancing algo")(
                "LoadBalancing.max_iterations",
                po::value<short>()->default_value(50),
                "Maximum number of iterations for load balancing algo")(
                "LoadBalancing.modulo", po::value<short>()->default_value(1),
                "Modulos or parameter to control how often clusters are "
                "recomputed during md")("LoadBalancing.output_file",
                po::value<string>()->default_value(""),
                "Output file for dumping cluster information in vtk format");

            // Hidden options, will be allowed in config file, but will not be
            // shown to the user.
            po::options_description hidden("Hidden options");
            hidden.add_options()("Quench.interval_print_residual",
                po::value<short>()->default_value(0),
                "Print interval for residual in quench (0 for none)")(
                "Quench.required_tol", po::value<float>()->default_value(1000.),
                "required tolerance (code stops if not reached)")(
                "Quench.compute_cond_Gram",
                po::value<bool>()->default_value(false),
                "Compute condition number of S during quench")(
                "Quench.min_Gram_eigenvalue",
                po::value<float>()->default_value(0.),
                "min. eigenvalue for Gram matrix")("Quench.step_length",
                po::value<float>()->default_value(-1.),
                "Step length for corrections in quench")("Quench.ortho_freq",
                po::value<short>()->default_value(1000),
                "Orthonormalization frequency")(
                "Quench.pair_mlwf_distance_threshold",
                po::value<float>()->default_value(4.),
                "Max. distance between pairs for MLWF transform")(
                "AOMM.kernel_radius", po::value<float>()->default_value(-1.),
                "Radius of kernel functions in AOMM algorithm")(
                "AOMM.threshold_factor", po::value<float>()->default_value(-1.),
                "Multiplicative factor of kernel radius to set threshold in "
                "AOMM projectors dropping")("Orbitals.type",
                po::value<string>()->default_value("NO"),
                "orbital type")("Orbitals.dotProduct",
                po::value<string>()->default_value("diagonal"),
                "orbital dot product type")("Orbitals.overallocate_factor",
                po::value<float>()->default_value(1.2),
                "safety factor to use for static allocation of orbitals")(
                "Potentials.filterPseudo",
                po::value<char>()->default_value('f'),
                "filter")("Poisson.solver",
                po::value<string>()->default_value("CG"), "solver")(
                "Poisson.rho0", po::value<float>()->default_value(0.0004),
                "continuum solvent: rho0")("Poisson.beta",
                po::value<float>()->default_value(1.3),
                "continuum solvent: beta")("Poisson.FDtype",
                po::value<string>()->default_value("Mehrstellen"), "FDtype")(
                "Poisson.nu1", po::value<short>()->default_value(2), "nu_1")(
                "Poisson.nu2", po::value<short>()->default_value(2), "nu_2")(
                "Poisson.max_steps", po::value<short>()->default_value(20),
                "max. nb. steps Poisson solver")("Poisson.max_steps_initial",
                po::value<short>()->default_value(20),
                "max. nb. steps Poisson solver in first solve")(
                "Poisson.max_levels", po::value<short>()->default_value(10),
                "max. nb. MG levels Poisson solver")("Poisson.reset",
                po::value<bool>()->default_value(false),
                "reset Hartree potential at each MD step")("ABPG.m",
                po::value<short>()->default_value(1),
                "History length for Anderson extrapolation")("ABPG.beta",
                po::value<float>()->default_value(1.),
                "beta for Anderson extrapolation")("NLCG.parallel_transport",
                po::value<bool>()->default_value(true),
                "Turn ON/OFF parallel transport algorithm")(
                "MD.extrapolation_type", po::value<short>()->default_value(1),
                "MD extrapolation type")("MD.compute_cond_Gram",
                po::value<bool>()->default_value(false),
                "Compute condition number of S at end of quench")(
                "MD.min_Gram_eigenvalue", po::value<float>()->default_value(0.),
                "min. eigenvalue for Gram matrix")(
                "ShortSightedInverse.spread_factor",
                po::value<float>()->default_value(2.),
                "Shortsighted spread factor")("ShortSightedInverse.tol",
                po::value<float>()->default_value(1.e-10),
                "Shortsighted tolerance")("ShortSightedInverse.krylov_dim",
                po::value<short>()->default_value(10),
                "Shortsighted Krylov space max. dimension")(
                "ShortSightedInverse.max_iterations",
                po::value<short>()->default_value(10),
                "Shortsighted max. number of GMRES iterations")(
                "ShortSightedInverse.ilu_type",
                po::value<string>()->default_value("ILU"),
                "Shortsighted ILU type: ILU or ILUT")(
                "ShortSightedInverse.ilu_drop_tol",
                po::value<float>()->default_value(1.e-5),
                "Shortsighted ILU dropping tolerance")(
                "ShortSightedInverse.ilu_filling_level",
                po::value<short>()->default_value(0),
                "Shortsighted ILU filling level")(
                "ShortSightedInverse.ilut_max_fill",
                po::value<int>()->default_value(10000),
                "Shortsighted max. filling for ILUT")("Coloring.algo",
                po::value<string>()->default_value("RLF"),
                "Coloring algorithm: RLF or Greedy")("Coloring.scope",
                po::value<string>()->default_value("local"),
                "Coloring scope: local or global")(
                "LocalizationRegions.min_distance",
                po::value<float>()->default_value(0.),
                "min. distance between Localization centers")(
                "LocalizationRegions.extrapolation_scheme",
                po::value<string>()->default_value("linear"),
                "Extrapolation order for localization centers")(
                "LocalizationRegions.computation",
                po::value<short>()->default_value(0),
                "Flag for computing new centers from extrapolated orbitals.")(
                "DensityMatrix.mixing", po::value<float>()->default_value(1.),
                "Mixing coefficient for Density Matrix")("DensityMatrix.solver",
                po::value<string>()->default_value("Mixing"),
                "Algorithm for updating Density Matrix: Mixing, MVP, HMVP")(
                "DensityMatrix.nb_inner_it",
                po::value<short>()->default_value(3),
                "Max. number of inner iterations in DM optimization")(
                "DensityMatrix.algo",
                po::value<string>()->default_value("Diagonalization"),
                "Algorithm for computing Density Matrix. "
                "Diagonalization or SP2.")("DensityMatrix.use_old",
                po::value<bool>()->default_value(true),
                "Start DM optimization with matrix of previous WF step")(
                "DensityMatrix.tol", po::value<float>()->default_value(1.e-7),
                "tolerance, used in iterative DM computation convergence "
                "criteria");

            po::options_description cmdline_options;
            cmdline_options.add(generic);

            po::options_description config_file_options;
            config_file_options.add(config).add(hidden);

            po::options_description visible("Allowed options");
            visible.add(generic).add(config);

            po::positional_options_description pd;
            pd.add("atomicCoordinates", -1);

            store(po::command_line_parser(argc, argv)
                      .options(cmdline_options)
                      .positional(pd)
                      .run(),
                vm);
            notify(vm);

            ifstream ifs(config_file.c_str());
            if (!ifs)
            {
                cout << "can not open config file: " << config_file << "\n";
                return 0;
            }
            else
            {
                store(parse_config_file(ifs, config_file_options), vm);
                notify(vm);
            }

            if (vm.count("help"))
            {
                if (onpe0) cout << visible << "\n";
                return 0;
            }

            if (vm.count("version"))
            {
                if (onpe0)
                {
#ifdef GITHASH
#define xstr(x) #x
#define LOG(x) cout << " MGmol: git_hash " << xstr(x) << endl;
                    LOG(GITHASH);
                    cout << endl;
#endif
                }
                return 0;
            }

            if (vm.count("check"))
            {
                tcheck = true;
            }
            if (vm.count("spin"))
            {
                total_spin = vm["spin"].as<float>();
                with_spin  = true;
                if (onpe0) cout << "Spin was set to " << total_spin << endl;
            }
            if (vm.count("atomicCoordinates"))
            {
                if (onpe0)
                    cout << "Input files is: "
                         << vm["atomicCoordinates"].as<vector<string>>()[0]
                         << "\n";
                input_file = vm["atomicCoordinates"].as<vector<string>>()[0];
            }
            else
            {
                if (vm["Restart.input_level"].as<short>() == 0)
                {
                    throw std::runtime_error(
                        "ERROR: No coordinates file provided!!!");
                }
            }

            if (onpe0) cout << "Spin: " << total_spin << "\n";

        } // try
        catch (exception& e)
        {
            cerr << e.what() << "\n";
            return 1;
        }

    MGmol_MPI::setup(MPI_COMM_WORLD, std::cout, with_spin);
    MGmol_MPI& mmpi      = *(MGmol_MPI::instance());
    MPI_Comm global_comm = mmpi.commGlobal();

    Control::setup(global_comm, with_spin, total_spin);
    Control& ct = *(Control::instance());

    ct.setOptions(vm);
    ct.sync();

    int ret = ct.checkOptions();
    if (ret < 0) return ret;

    mmpi.bcastGlobal(input_file);
    mmpi.bcastGlobal(lrs_filename);

#ifdef _OPENMP
    if (onpe0)
    {
        cout << " " << omp_get_max_threads() << " thread"
             << (omp_get_max_threads() > 1 ? "s " : " ");
        cout << "active" << endl << endl;
    }
    omp_set_nested(0);
    if (omp_get_nested())
    {
        cerr << "Nested parallelism not allowed" << endl;
        return 1;
    }
#endif

    // Enter main scope
    {
        // setup standard output and error
        sout = &std::cout;
        serr = &std::cerr;

        MGmolInterface* mgmol;
        if (ct.isLocMode())
            mgmol = new MGmol<LocGridOrbitals>(global_comm, *MPIdata::sout);
        else
            mgmol
                = new MGmol<ExtendedGridOrbitals>(global_comm, *MPIdata::sout);

        unsigned ngpts[3]    = { ct.ngpts_[0], ct.ngpts_[1], ct.ngpts_[2] };
        double origin[3]     = { ct.ox_, ct.oy_, ct.oz_ };
        const double cell[3] = { ct.lx_, ct.ly_, ct.lz_ };
        Mesh::setup(mmpi.commSpin(), ngpts, origin, cell, ct.lap_type);

        mgmol->setupFromInput(input_file);
        if (ct.restart_info < 3 || !ct.isLocMode())
            mgmol->setupLRsFromInput(lrs_filename);
        mgmol->setupConstraintsFromInput(constraints_filename);

        ct.checkNLrange();

        LocGridOrbitals::setDotProduct(ct.dot_product_type);

        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();

        if (!ct.short_sighted)
        {
            MatricesBlacsContext::instance().setup(mmpi.commSpin(), ct.numst);

            dist_matrix::DistMatrix<DISTMATDTYPE>::setBlockSize(64);

            dist_matrix::DistMatrix<DISTMATDTYPE>::setDefaultBlacsContext(
                MatricesBlacsContext::instance().bcxt());

            ReplicatedWorkSpace<double>::instance().setup(ct.numst);

            dist_matrix::SparseDistMatrix<
                DISTMATDTYPE>::setNumTasksPerPartitioning(128);

            MGmol_MPI& mmpi = *(MGmol_MPI::instance());
            int npes        = mmpi.size();
            setSparseDistMatriConsolidationNumber(npes);
        }

#ifdef USE_DIS_MAT
        if (myPEenv.color() > 0)
        {
            cerr << "Code should be called with " << myPEenv.n_mpi_tasks()
                 << " MPI tasks only" << endl;
            ct.global_exit(2);
        }

#endif

        assert(myPEenv.color() == 0);

        if (myPEenv.color() == 0)
        {
            assert(ct.getMGlevels() >= -1);
            if (ct.getMGlevels() >= 0 && ct.getPrecondType() % 10 == 0)
            {
                const pb::Grid& mygrid = mymesh->grid();

                if ((mygrid.dim(0) % (1 << ct.getMGlevels())) != 0)
                {
                    cerr << "main: mygrid.dim(0)=" << mygrid.dim(0)
                         << " not evenly divisible by "
                         << "2^(ct.getMGlevels())=" << (1 << ct.getMGlevels())
                         << endl;
                    return -1;
                }
                if ((mygrid.dim(1) % (1 << ct.getMGlevels())) != 0)
                {
                    cerr << "main: mygrid.dim(1)=" << mygrid.dim(1)
                         << " not evenly divisible by "
                         << "2^(ct.getMGlevels())=" << (1 << ct.getMGlevels())
                         << endl;
                    return -1;
                }
                if ((mygrid.dim(2) % (1 << ct.getMGlevels())) != 0)
                {
                    cerr << "main: mygrid.dim(2)=" << mygrid.dim(2)
                         << " not evenly divisible by "
                         << "2^(ct.getMGlevels())=" << (1 << ct.getMGlevels())
                         << endl;
                    return -1;
                }
            }

            myPEenv.barrier(); // wait to see if everybody is OK before
                               // continuing...

            if (!tcheck)
            {
#ifdef DEBUG
                *MPIdata::sout << " Run begins on processor "
                               << myPEenv.mytask() << endl;
#endif

#ifdef USE_CNR
                /* use conditional numerical reproducibility for MKL */
                int my_cbwr_branch = mkl_cbwr_get_auto_branch();
                if (mkl_cbwr_set(my_cbwr_branch) != MKL_CBWR_SUCCESS)
                {
                    printf("Error in setting MKL_CBWR_BRANCH! Aborting…\n");
                    return (-1);
                }
#endif
                mgmol->run();
            }
            else
            {
                *MPIdata::sout << " Input parameters OK\n";
            }
        }
        delete mgmol;

        if (!ct.short_sighted)
        {
            // need to destroy any MPI based object befor calling MPI_Finalize
            MatricesBlacsContext::instance().clear();
        }
    } // close main scope

    // release memory for static arrays
    PackedCommunicationBuffer::deleteStorage();
    Mesh::deleteInstance();
    Control::deleteInstance();
    MGmol_MPI::deleteInstance();

#ifdef HAVE_MAGMA
    // Delete the data in the singleton before finalizing magma
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();
    magma_singleton.free();

    magmalog = magma_finalize();

    if (magmalog == MAGMA_SUCCESS)
    {
        std::cout << "MAGMA Finalize: success" << std::endl;
    }
    else
    {
        if (magmalog == MAGMA_ERR_UNKNOWN)
            std::cout << "MAGMA Finalize: unknown error" << std::endl;
        if (magmalog == MAGMA_ERR_HOST_ALLOC)
            std::cout << "MAGMA FINALIZE: fails host alloc" << std::endl;
        return 1;
    }

#endif

    mpirc = MPI_Finalize();
    if (mpirc != MPI_SUCCESS)
    {
        cerr << "MPI Finalize failed!!!" << endl;
    }

    time_t tt;
    time(&tt);
    if (onpe0) cout << " Run ended at " << ctime(&tt) << endl;

    // MemTrack::TrackDumpBlocks();

    //    MemTrack::TrackListMemoryUsage();

    return 0;
}
