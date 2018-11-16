// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "ReplicatedWorkSpace.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

// pb
#include "MGmol_blas1.h"
#include "mputils.h"

#include "Control.h"
#include "MGmol_MPI.h"

Timer ReplicatedWorkSpace::mpisum_tm_("ReplicatedWorkSpace::mpisum");

// allows zero or one object
ReplicatedWorkSpace& ReplicatedWorkSpace::instance()
{
    static ReplicatedWorkSpace instance;
    return instance;
}

void ReplicatedWorkSpace::mpiBcastSquareMatrix()
{
#ifdef USE_MPI
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int n2          = ndim_ * ndim_;
    mmpi.bcast(square_matrix_, n2, 0);
//    MPI_Bcast(square_matrix_, n2, MPI_DOUBLE, 0, mmpi.commSpin());
#endif
    return;
}

void ReplicatedWorkSpace::setUpperTriangularSquareMatrixToZero()
{
    for (int i = 0; i < ndim_; i++)
    {
        for (int j = i + 1; j < ndim_; j++)
        {
            square_matrix_[i + ndim_ * j] = 0.;
        }
    }
}
