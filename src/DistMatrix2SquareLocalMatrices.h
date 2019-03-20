// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DistMatrix2SquareLocalMatrices_H
#define MGMOL_DistMatrix2SquareLocalMatrices_H

#include "DistMatrix.h"
#include "SubMatricesIndexing.h"
#include "SubMatrices.h"
#include "SquareLocalMatrices.h"
#include "Timer.h"

#include <memory>
#include <vector>

#include <mpi.h>

class DistMatrix2SquareLocalMatrices
{
    static Timer convert_tm_;

    std::unique_ptr< dist_matrix::SubMatricesIndexing<DISTMATDTYPE> > submat_indexing_;

    std::unique_ptr< dist_matrix::SubMatrices<DISTMATDTYPE> > submatWork_;

public:

    DistMatrix2SquareLocalMatrices(
        MPI_Comm comm,
        const std::vector<std::vector<int>>& gids,
        const dist_matrix::DistMatrix<DISTMATDTYPE>& empty_mat)
        : submat_indexing_( new dist_matrix::SubMatricesIndexing<DISTMATDTYPE>(
                            gids, comm, empty_mat) ),
          submatWork_( new dist_matrix::SubMatrices<DISTMATDTYPE>("Work",
                       gids, comm, empty_mat, *submat_indexing_) )
    {
    }

    ~DistMatrix2SquareLocalMatrices()
    {
    }

    void convert(const dist_matrix::DistMatrix<DISTMATDTYPE>& dmat,
                SquareLocalMatrices<MATDTYPE>& lmat);

    static void printTimers(std::ostream& os)
    {
        convert_tm_.print(os);
    }
};

#endif

