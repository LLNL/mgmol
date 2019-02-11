// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ReplicatedWorkSpace.h"

#include "MGmol_blas1.h"
#include "MGmol_MPI.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

template <class T>
Timer ReplicatedWorkSpace<T>::mpisum_tm_("ReplicatedWorkSpace::mpisum");

// allows zero or one object
template <class T>
ReplicatedWorkSpace<T>& ReplicatedWorkSpace<T>::instance()
{
    static ReplicatedWorkSpace<T> instance;
    return instance;
}

template <class T>
void ReplicatedWorkSpace<T>::mpiBcastSquareMatrix()
{
#ifdef USE_MPI
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int n2          = ndim_ * ndim_;
    mmpi.bcast(square_matrix_, n2, 0);
#endif
    return;
}

template <class T>
void ReplicatedWorkSpace<T>::setUpperTriangularSquareMatrixToZero()
{
    for (int i = 0; i < ndim_; i++)
    {
        for (int j = i + 1; j < ndim_; j++)
        {
            square_matrix_[i + ndim_ * j] = 0.;
        }
    }
}

template class ReplicatedWorkSpace<double>;

