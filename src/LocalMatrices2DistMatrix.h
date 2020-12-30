// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_LocalMatrices2DistMatrix_H
#define MGMOL_LocalMatrices2DistMatrix_H

#include "LocalMatrices.h"
#include "SparseDistMatrix.h"

#include <vector>

// Add matrix elements corresponding to subdomains at their right place
// into distributed matrix
// Important Note: Neglect contributions smaller than tol!
// (may lead to results dependent on number of CPUs)

class LocalMatrices2DistMatrix
{
private:
    static LocalMatrices2DistMatrix* pinstance_;

    static MPI_Comm comm_;

    static std::vector<std::vector<int>> global_indexes_;

    static double tol_mat_elements;

    static short sparse_distmatrix_nb_tasks_per_partitions_;

public:
    static LocalMatrices2DistMatrix* instance()
    {
        if (pinstance_ == nullptr)
        {
            pinstance_ = new LocalMatrices2DistMatrix();
        }
        return pinstance_;
    }

    LocalMatrices2DistMatrix() { assert(comm_ != MPI_COMM_NULL); }

    static void setup(MPI_Comm comm, const std::vector<std::vector<int>>& gids)
    {
        comm_           = comm;
        global_indexes_ = gids;
    }

    template <typename ScalarType>
    void convert(const LocalMatrices<ScalarType, MemorySpace::Host>& lmat,
        dist_matrix::SparseDistMatrix<ScalarType>& dst, const int numst,
        const double tol = tol_mat_elements) const;

    template <class ScalarType>
    void accumulate(const LocalMatrices<ScalarType, MemorySpace::Host>& lmat,
        dist_matrix::DistMatrix<ScalarType>& dst,
        const double tol = tol_mat_elements) const;
};

#endif
