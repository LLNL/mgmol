// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ReplicatedWorkSpace.h"

#include "MGmol_MPI.h"
#include "MGmol_blas1.h"

#include <cassert>
#include <mpi.h>

template <class ScalarType>
Timer ReplicatedWorkSpace<ScalarType>::mpisum_tm_(
    "ReplicatedWorkSpace::mpisum");

// allows zero or one object
template <class ScalarType>
ReplicatedWorkSpace<ScalarType>& ReplicatedWorkSpace<ScalarType>::instance()
{
    static ReplicatedWorkSpace<ScalarType> instance;
    return instance;
}

template <class ScalarType>
void ReplicatedWorkSpace<ScalarType>::mpiBcastSquareMatrix()
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int n2          = ndim_ * ndim_;
    mmpi.bcast(square_matrix_, n2, 0);
    return;
}

template <class ScalarType>
void ReplicatedWorkSpace<ScalarType>::setUpperTriangularSquareMatrixToZero()
{
    for (int i = 0; i < ndim_; i++)
    {
        for (int j = i + 1; j < ndim_; j++)
        {
            square_matrix_[i + ndim_ * j] = 0.;
        }
    }
}

#ifdef MGMOL_USE_SCALAPACK
template <class ScalarType>
void ReplicatedWorkSpace<ScalarType>::initSquareMatrix(
    const dist_matrix::DistMatrix<ScalarType>& distmat)
{
    distmat.allgather(square_matrix_, ndim_);
}
#endif

template <class ScalarType>
void ReplicatedWorkSpace<ScalarType>::initSquareMatrix(
    const ReplicatedMatrix& mat)
{
    assert(square_matrix_ != nullptr);
    assert(ndim_ > 0);

    mat.get(square_matrix_, ndim_);
}

template class ReplicatedWorkSpace<double>;
