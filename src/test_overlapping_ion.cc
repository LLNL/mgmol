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
#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "mgmol_run.h"

#include <cassert>
#include <iostream>
#include <time.h>
#include <vector>
#include <random>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

template <class OrbitalsType>
void testOverlappingIons(MGmolInterface *mgmol_)
{
    /* random number generator */
    static std::random_device rd;  // Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd(){}
    static std::uniform_real_distribution<> dis(0.0, 1.0);

    MGmol<OrbitalsType> *mgmol = static_cast<MGmol<OrbitalsType> *>(mgmol_);
    std::shared_ptr<Ions> ions = mgmol->getIons();

    /* 3 fictitious ion configurations */
    const int num_snap = 3;
    const std::vector<Ion*> &local_ions(ions->local_ions());
    const int nlocal_ions = local_ions.size();
    std::vector<std::vector<std::vector<double>>> cfgs(num_snap);
    for (int idx = 0; idx < num_snap; idx++)
    {
        cfgs[idx].resize(nlocal_ions);

        for (int k = 0; k < nlocal_ions; k++)
        {
            cfgs[idx][k].resize(3);

            for (int d = 0; d < 3; d++)
            {
                double orig_position = local_ions[k]->position(d);
                /* hope this does not go beyond the domain.. */
                cfgs[idx][k][d] = orig_position * (0.9 + 0.2 * dis(gen));
            }
        }
    }

    /* Save overlappingVL ions from each fictitious ion configuration */
    std::vector<std::vector<std::vector<double>>> fom_overlapping_ions(num_snap);
    for (int idx = 0; idx < num_snap; idx++)
    {
        /* set ion positions */
        for (int i = 0; i < nlocal_ions; i++)
            local_ions[i]->setPosition(cfgs[idx][i][0], cfgs[idx][i][1], cfgs[idx][i][2]);
        ions->setup();

        fom_overlapping_ions[idx].resize(ions->overlappingNL_ions().size());
        for (int k = 0; k < ions->overlappingNL_ions().size(); k++)
        {
            fom_overlapping_ions[idx][k].resize(3);

            for (int d = 0; d < 3; d++)
                fom_overlapping_ions[idx][k][d] = ions->overlappingNL_ions()[k]->position(d);
        }
    }

    const int test_idx = 1;
    /* set ion positions */
    for (int i = 0; i < nlocal_ions; i++)
        local_ions[i]->setPosition(cfgs[test_idx][i][0], cfgs[test_idx][i][1], cfgs[test_idx][i][2]);
    ions->setup();

    if (fom_overlapping_ions[test_idx].size() != ions->overlappingNL_ions().size())
    {
        std::cerr << "overlapping ion number is different!!!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    for (int k = 0; k < ions->overlappingNL_ions().size(); k++)
    {
        for (int d = 0; d < 3; d++)
        {
            if (abs(fom_overlapping_ions[test_idx][k][d] - ions->overlappingNL_ions()[k]->position(d)) >= 1.0e-12)
            {
                std::cerr << "overlapping ion coordinates are different!!!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }
            
    }
}

int main(int argc, char** argv)
{
    int mpirc = MPI_Init(&argc, &argv);
    if (mpirc != MPI_SUCCESS)
    {
        std::cerr << "MPI Initialization failed!!!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    MPI_Comm comm = MPI_COMM_WORLD;

    /*
     * Initialize general things, like magma, openmp, IO, ...
     */
    mgmol_init(comm);

    /*
     * read runtime parameters
     */
    std::string input_filename("");
    std::string lrs_filename;
    std::string constraints_filename("");

    float total_spin = 0.;
    bool with_spin   = false;

    po::variables_map vm;

    // read from PE0 only
    if (MPIdata::onpe0)
    {
        read_config(argc, argv, vm, input_filename, lrs_filename,
            constraints_filename, total_spin, with_spin);
    }

    MGmol_MPI::setup(comm, std::cout, with_spin);
    MGmol_MPI& mmpi      = *(MGmol_MPI::instance());
    MPI_Comm global_comm = mmpi.commGlobal();

    /*
     * Setup control struct with run time parameters
     */
    Control::setup(global_comm, with_spin, total_spin);
    Control& ct = *(Control::instance());

    ct.setOptions(vm);

    int ret = ct.checkOptions();
    if (ret < 0) return ret;

    mmpi.bcastGlobal(input_filename);
    mmpi.bcastGlobal(lrs_filename);

    // Enter main scope
    {
        MGmolInterface* mgmol;
        if (ct.isLocMode())
            mgmol = new MGmol<LocGridOrbitals>(global_comm, *MPIdata::sout,
                input_filename, lrs_filename, constraints_filename);
        else
            mgmol = new MGmol<ExtendedGridOrbitals>(global_comm, *MPIdata::sout,
                input_filename, lrs_filename, constraints_filename);

        mgmol->setup();

        if (ct.isLocMode())
            testOverlappingIons<LocGridOrbitals>(mgmol);
        else
            testOverlappingIons<ExtendedGridOrbitals>(mgmol);

        delete mgmol;

    } // close main scope

    mgmol_finalize();

    mpirc = MPI_Finalize();
    if (mpirc != MPI_SUCCESS)
    {
        std::cerr << "MPI Finalize failed!!!" << std::endl;
    }

    time_t tt;
    time(&tt);
    if (onpe0) std::cout << " Run ended at " << ctime(&tt) << std::endl;

    return 0;
}
