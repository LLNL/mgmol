// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SquareSubMatrix2DistMatrix.h"
#include "MGmol_MPI.h"

SquareSubMatrix2DistMatrix* SquareSubMatrix2DistMatrix::pinstance_ = nullptr;
double SquareSubMatrix2DistMatrix::tol_mat_elements                = 1.e-14;
short SquareSubMatrix2DistMatrix::sparse_distmatrix_nb_tasks_per_partitions_
    = 256;

template <class T>
void SquareSubMatrix2DistMatrix::convert(const SquareSubMatrix<T>& src,
    dist_matrix::SparseDistMatrix<T>& dst, const double tol) const
{
    dst.addData(src, tol);
}

// Sum up all the local contributions (in SquareSubMatrix) into
// one DistMatrix
template <class T>
void SquareSubMatrix2DistMatrix::accumulate(const SquareSubMatrix<T>& src,
    dist_matrix::DistMatrix<T>& dst, const double tol) const
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sm(
        comm, dst, sparse_distmatrix_nb_tasks_per_partitions_);

    // convert into a SparseDistMatrix
    convert(src, sm, tol);

    // accumulate into a DistMatrix
    sm.parallelSumToDistMatrix();
}

template void SquareSubMatrix2DistMatrix::convert(
    const SquareSubMatrix<DISTMATDTYPE>& src,
    dist_matrix::SparseDistMatrix<DISTMATDTYPE>& dst, const double tol) const;
template void SquareSubMatrix2DistMatrix::accumulate(
    const SquareSubMatrix<DISTMATDTYPE>& src,
    dist_matrix::DistMatrix<DISTMATDTYPE>& dst, const double tol) const;
