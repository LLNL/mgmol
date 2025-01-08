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
#include "Potentials.h"

#include <cassert>
#include <iostream>
#include <time.h>
#include <vector>
#include <random>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

template <class OrbitalsType>
void testRhoRestart(MGmolInterface *mgmol_)
{
    /* random number generator */
    static std::random_device rd;  // Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd(){}
    static std::uniform_real_distribution<> dis(0.0, 1.0);

    MGmol_MPI& mmpi      = *(MGmol_MPI::instance());
    const int rank = mmpi.mypeGlobal();
    const int nprocs = mmpi.size();

    Control& ct = *(Control::instance());
    MGmol<OrbitalsType> *mgmol = static_cast<MGmol<OrbitalsType> *>(mgmol_);
    std::shared_ptr<Rho<OrbitalsType>> rho = mgmol->getRho();

    /* load a restart file */
    mgmol->loadRestartFile(ct.restart_file);

    /* save density from the restart file to elsewhere */
    std::vector<RHODTYPE> rho0(rho->rho_[0].size());
    rho0 = rho->rho_[0];

    /* recompute rho from the orbital */
    rho->update(*mgmol->getOrbitals());

    /* check if the recomputed density is the same */
    for (int d = 0; d < rho0.size(); d++)
    {
        double error = abs(rho0[d] - rho->rho_[0][d]) / abs(rho0[d]);
        if (error > 1e-10)
        {
            printf("rank %d, rho[%d]=%.15e, rho0[%d]=%.15e\n", rank, d, rho->rho_[0][d], d, rho0[d]);
            std::cerr << "Density is inconsistent!!!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 0);
        }
    }
}

template <class OrbitalsType>
void testPotRestart(MGmolInterface *mgmol_)
{
    /* random number generator */
    static std::random_device rd;  // Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd(){}
    static std::uniform_real_distribution<> dis(0.0, 1.0);

    MGmol_MPI& mmpi      = *(MGmol_MPI::instance());
    const int rank = mmpi.mypeGlobal();
    const int nprocs = mmpi.size();

    Control& ct = *(Control::instance());
    MGmol<OrbitalsType> *mgmol = static_cast<MGmol<OrbitalsType> *>(mgmol_);
    Potentials& pot = mgmol->getHamiltonian()->potential();

    /* load a restart file */
    mgmol->loadRestartFile(ct.restart_file);

    /* save potential from the restart file to elsewhere */
    std::vector<POTDTYPE> vh0(pot.size());
    POTDTYPE *d_vhrho = pot.vh_rho();
    for (int d = 0; d < vh0.size(); d++)
        vh0[d] = d_vhrho[d];

    /* recompute potential */
    std::shared_ptr<Ions> ions = mgmol->getIons();
    mgmol->update_pot(*ions);

    /* check if the recomputed potential is the same */
    d_vhrho = pot.vh_rho();
    for (int d = 0; d < vh0.size(); d++)
    {
        double error = abs(vh0[d] - d_vhrho[d]) / abs(vh0[d]);
        if (error > 1e-10)
        {
            printf("rank %d, vh_rho[%d]=%.15e, vh_rho0[%d]=%.15e\n", rank, d, d_vhrho[d], d, vh0[d]);
            std::cerr << "Potential is inconsistent!!!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 0);
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
        {
            testRhoRestart<LocGridOrbitals>(mgmol);
            testPotRestart<LocGridOrbitals>(mgmol);
        }
        else
        {
            testRhoRestart<ExtendedGridOrbitals>(mgmol);
            testPotRestart<ExtendedGridOrbitals>(mgmol);
        }

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
