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

    MGmol_MPI& mmpi      = *(MGmol_MPI::instance());
    const int rank = mmpi.mypeGlobal();
    const int nprocs = mmpi.size();

    MGmol<OrbitalsType> *mgmol = static_cast<MGmol<OrbitalsType> *>(mgmol_);
    std::shared_ptr<Ions> ions = mgmol->getIons();

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    /* get the extent of global domain */
    const double origin[3] = { mygrid.origin(0), mygrid.origin(1), mygrid.origin(2) };
    const double lattice[3] = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    if (rank == 0)
    {
        printf("origin: (%.3e, %.3e, %.3e)\n", origin[0], origin[1], origin[2]);
        printf("lattice: (%.3e, %.3e, %.3e)\n", lattice[0], lattice[1], lattice[2]);
    }

    /* get global atomic numbers */
    const int num_ions = ions->getNumIons();
    std::vector<short> atnumbers(num_ions);
    ions->getAtomicNumbers(atnumbers);

    /* 3 fictitious ion configurations */
    const int num_snap = 3;
    const std::vector<Ion*> &local_ions(ions->local_ions());
    const int nlocal_ions = local_ions.size();
    std::vector<std::vector<double>> cfgs(num_snap);
    for (int idx = 0; idx < num_snap; idx++)
    {
        cfgs[idx].resize(3 * num_ions);
        if (rank == 0)
            for (int k = 0; k < num_ions; k++)
                for (int d = 0; d < 3; d++)
                    cfgs[idx][3 * k + d] = origin[d] + lattice[d] * dis(gen);
        
        mmpi.bcastGlobal(cfgs[idx].data(), 3 * num_ions, 0);
    }

    /* Save overlappingVL ions from each fictitious ion configuration */
    std::vector<std::vector<std::vector<double>>> fom_overlapping_ions(num_snap);
    for (int idx = 0; idx < num_snap; idx++)
    {
        /* set ion positions */
        ions->setPositions(cfgs[idx], atnumbers);

        fom_overlapping_ions[idx].resize(ions->overlappingNL_ions().size());
        for (int k = 0; k < ions->overlappingNL_ions().size(); k++)
        {
            fom_overlapping_ions[idx][k].resize(3);

            for (int d = 0; d < 3; d++)
                fom_overlapping_ions[idx][k][d] = ions->overlappingNL_ions()[k]->position(d);
        }
    }

    std::uniform_int_distribution<> distrib(0, num_snap-1);
    int test_idx = distrib(gen);
    mmpi.bcastGlobal(&test_idx);
    if (rank == 0) printf("test index: %d\n", test_idx);

    /* set ion positions */
    ions->setPositions(cfgs[test_idx], atnumbers);

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
                printf("rank %d, local overlap ion %d: (%.3e, %.3e, %.3e) =?= (%.3e, %.3e, %.3e)\n",
                        rank, k, fom_overlapping_ions[test_idx][k][0], fom_overlapping_ions[test_idx][k][1], fom_overlapping_ions[test_idx][k][2],
                        ions->overlappingNL_ions()[k]->position(0), ions->overlappingNL_ions()[k]->position(1), ions->overlappingNL_ions()[k]->position(2));
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
