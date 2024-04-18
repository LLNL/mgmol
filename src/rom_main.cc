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

    // release memory for static arrays
    // PackedCommunicationBuffer::deleteStorage();
    // Mesh::deleteInstance();
    Control::deleteInstance();
    // MGmol_MPI::deleteInstance();

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
