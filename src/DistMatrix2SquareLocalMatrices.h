// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DistMatrix2SquareLocalMatrices_H
#define MGMOL_DistMatrix2SquareLocalMatrices_H

#include "DistMatrix.h"
#include "SquareLocalMatrices.h"
#include "SubMatrices.h"
#include "SubMatricesIndexing.h"
#include "Timer.h"

#include <memory>
#include <vector>

#include <mpi.h>

class DistMatrix2SquareLocalMatrices
{
    static DistMatrix2SquareLocalMatrices* pinstance_;

    static Timer convert_tm_;

    static MPI_Comm comm_;

    static std::unique_ptr<dist_matrix::SubMatricesIndexing<DISTMATDTYPE>>
        submat_indexing_;
    static std::unique_ptr<dist_matrix::SubMatrices<DISTMATDTYPE>> submatWork_;

public:
    static DistMatrix2SquareLocalMatrices* instance()
    {
        if (pinstance_ == nullptr)
        {
            pinstance_ = new DistMatrix2SquareLocalMatrices();
        }
        return pinstance_;
    }

    DistMatrix2SquareLocalMatrices()
    {
        assert(comm_ != MPI_COMM_NULL);
        assert(submat_indexing_);
        assert(submatWork_);
    }

    static void setup(MPI_Comm comm, const std::vector<std::vector<int>>& gids,
        const dist_matrix::DistMatrix<DISTMATDTYPE>& empty_mat)
    {
        comm_ = comm;
        submat_indexing_.reset(
            new dist_matrix::SubMatricesIndexing<DISTMATDTYPE>(
                gids, comm, empty_mat));
        submatWork_.reset(new dist_matrix::SubMatrices<DISTMATDTYPE>(
            "Work", gids, comm, empty_mat, *submat_indexing_));
    }

    ~DistMatrix2SquareLocalMatrices() { }

    void convert(const dist_matrix::DistMatrix<DISTMATDTYPE>& dmat,
        SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& lmat);
    const dist_matrix::SubMatrices<DISTMATDTYPE>& convert(
        const dist_matrix::DistMatrix<DISTMATDTYPE>& dmat);

    static void printTimers(std::ostream& os) { convert_tm_.print(os); }
};

#endif
