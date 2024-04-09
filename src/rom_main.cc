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

#include "librom.h"

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

    CAROM::Vector librom_vector(10, false);

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
