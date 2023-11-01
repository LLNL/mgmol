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
#include <cassert>
#include <iostream>
#include <iterator>
#include <vector>

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
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " "));
    return os;
}

int main(int argc, char** argv)
{
    // change handling of memory allocation errors
    std::set_new_handler(noMoreMemory);

    std::cout.sync_with_stdio();

    int mpirc = MPI_Init(&argc, &argv);
    if (mpirc != MPI_SUCCESS)
    {
        std::cerr << "MPI Initialization failed!!!" << std::endl;
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

    std::string input_file("");
    std::string lrs_filename;
    std::string constraints_filename("");
    bool tcheck = false;

    float total_spin = 0.;
    bool with_spin   = false;

    po::variables_map vm;

    // use configure file if it can be found
    // std::string config_filename("mgmol.cfg");

    // read options from PE0 only
    if (onpe0)
    {
        read_config(argc, argv, vm, input_file, lrs_filename,
            constraints_filename, total_spin, with_spin, tcheck);
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
        std::cout << " " << omp_get_max_threads() << " thread"
                  << (omp_get_max_threads() > 1 ? "s " : " ");
        std::cout << "active" << std::endl << std::endl;
    }
    omp_set_nested(0);
    if (omp_get_nested())
    {
        std::cerr << "Nested parallelism not allowed" << std::endl;
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
            std::cerr << "Code should be called with " << myPEenv.n_mpi_tasks()
                      << " MPI tasks only" << std::endl;
            ct.global_exit(2);
        }

        assert(myPEenv.color() == 0);

        {
            assert(ct.getMGlevels() >= -1);
            if (ct.withPreconditioner())
            {
                const pb::Grid& mygrid = mymesh->grid();

                if ((mygrid.dim(0) % (1 << ct.getMGlevels())) != 0)
                {
                    std::cerr
                        << "main: mygrid.dim(0)=" << mygrid.dim(0)
                        << " not evenly divisible by "
                        << "2^(ct.getMGlevels())=" << (1 << ct.getMGlevels())
                        << std::endl;
                    return -1;
                }
                if ((mygrid.dim(1) % (1 << ct.getMGlevels())) != 0)
                {
                    std::cerr
                        << "main: mygrid.dim(1)=" << mygrid.dim(1)
                        << " not evenly divisible by "
                        << "2^(ct.getMGlevels())=" << (1 << ct.getMGlevels())
                        << std::endl;
                    return -1;
                }
                if ((mygrid.dim(2) % (1 << ct.getMGlevels())) != 0)
                {
                    std::cerr
                        << "main: mygrid.dim(2)=" << mygrid.dim(2)
                        << " not evenly divisible by "
                        << "2^(ct.getMGlevels())=" << (1 << ct.getMGlevels())
                        << std::endl;
                    return -1;
                }
            }

            myPEenv.barrier(); // wait to see if everybody is OK before
                               // continuing...

            if (!tcheck)
            {
#ifdef DEBUG
                *MPIdata::sout << " Run begins on processor "
                               << myPEenv.mytask() << std::endl;
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
        std::cerr << "MPI Finalize failed!!!" << std::endl;
    }

    time_t tt;
    time(&tt);
    if (onpe0) std::cout << " Run ended at " << ctime(&tt) << std::endl;

    // MemTrack::TrackDumpBlocks();

    //    MemTrack::TrackListMemoryUsage();

    return 0;
}
