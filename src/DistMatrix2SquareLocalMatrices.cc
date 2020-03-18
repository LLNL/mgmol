// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DistMatrix2SquareLocalMatrices.h"

DistMatrix2SquareLocalMatrices* DistMatrix2SquareLocalMatrices::pinstance_
    = nullptr;
MPI_Comm DistMatrix2SquareLocalMatrices::comm_ = MPI_COMM_NULL;
std::unique_ptr<dist_matrix::SubMatricesIndexing<DISTMATDTYPE>>
    DistMatrix2SquareLocalMatrices::submat_indexing_;
std::unique_ptr<dist_matrix::SubMatrices<DISTMATDTYPE>>
    DistMatrix2SquareLocalMatrices::submatWork_;

Timer DistMatrix2SquareLocalMatrices::convert_tm_(
    "DistMatrix2SquareLocalMatrices::convert");

void DistMatrix2SquareLocalMatrices::convert(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& dmat,
    SquareLocalMatrices<MATDTYPE>& lmat)
{
    convert_tm_.start();

    submatWork_->gather(dmat);

    const short nd = lmat.subdiv();
    const int dim  = lmat.n();

    for (short i = 0; i < nd; i++)
    {
        DISTMATDTYPE* local = lmat.getSubMatrix(i);

        for (int ii = 0; ii < dim; ii++)
            for (int jj = 0; jj < dim; jj++)
                local[ii + dim * jj] = submatWork_->val(ii, jj, i);
    }

    convert_tm_.stop();
}

const dist_matrix::SubMatrices<DISTMATDTYPE>&
DistMatrix2SquareLocalMatrices::convert(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& dmat)
{
    submatWork_->gather(dmat);

    return *submatWork_;
}
