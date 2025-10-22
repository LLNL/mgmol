// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ReplicatedMatrix2SquareLocalMatrices.h"

ReplicatedMatrix2SquareLocalMatrices*
    ReplicatedMatrix2SquareLocalMatrices::pinstance_
    = nullptr;
std::vector<std::vector<int>>
    ReplicatedMatrix2SquareLocalMatrices::global_indexes_;

Timer ReplicatedMatrix2SquareLocalMatrices::convert_tm_(
    "ReplicatedMatrix2SquareLocalMatrices::convert");

void ReplicatedMatrix2SquareLocalMatrices::convert(const ReplicatedMatrix& rmat,
    SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& lmat)
{
    convert_tm_.start();

    const short nd = lmat.nmat();
    const int dim  = lmat.n();
    const int nst  = rmat.ld();

    for (short i = 0; i < nd; i++)
    {
        double* dst = lmat.getSubMatrix(i);
        double* src = rmat.data();
        for (int jj = 0; jj < dim; jj++)
        {
            const int st2 = global_indexes_[i][jj];
            if (st2 != -1)
            {
                for (int ii = 0; ii < dim; ii++)
                {
                    const int st1 = global_indexes_[i][ii];
                    if (st1 != -1)
                    {
                        dst[ii + dim * jj] = src[st1 + nst * st2];
                    }
                }
            }
        }
    }

    convert_tm_.stop();
}
