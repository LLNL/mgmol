// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
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
#include "rom_workflows.h"

//#include "MemTrack.h"

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
    std::cout << "ROM does not support MAGMA!\n" << std::endl;
    return 1;
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
            setupGenericOption(generic, config_file, lrs_filename);

            // Declare a group of options (with default when appropriate) that
            // will be allowed in config file
            po::options_description config("Configuration");
            setupConfigOption(config, constraints_filename);

            // Hidden options, will be allowed in config file, but will not be
            // shown to the user.
            po::options_description hidden("Hidden options");
            setupHiddenOption(hidden);

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

        if (ct.isLocMode() || ct.init_loc == 1) mgmol->setupLRs(lrs_filename);

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
#ifdef HAVE_MAGMA
        ReplicatedMatrix::setMPIcomm(mmpi.commSpin());
#endif

        if (myPEenv.color() > 0)
        {
            cerr << "Code should be called with " << myPEenv.n_mpi_tasks()
                 << " MPI tasks only" << endl;
            ct.global_exit(2);
        }

        assert(myPEenv.color() == 0);

        if (myPEenv.color() == 0)
        {
            assert(ct.getMGlevels() >= -1);
            if (ct.withPreconditioner())
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

            // mgmol->run();
            mgmol->setup();
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
