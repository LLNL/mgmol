// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_SquareSubMatrix2DistMatrix_H
#define MGMOL_SquareSubMatrix2DistMatrix_H

#include "SparseDistMatrix.h"
#include "SquareSubMatrix.h"

#include <vector>

// Add matrix elements corresponding to subdomains at their right place
// into distributed matrix
// Important Note: Neglect contributions smaller than tol!
// (may lead to results dependent on number of CPUs)

class SquareSubMatrix2DistMatrix
{
private:
    static SquareSubMatrix2DistMatrix* pinstance_;

    static double tol_mat_elements;

    static short sparse_distmatrix_nb_tasks_per_partitions_;

public:
    static SquareSubMatrix2DistMatrix* instance()
    {
        if (pinstance_ == nullptr)
        {
            pinstance_ = new SquareSubMatrix2DistMatrix();
        }
        return pinstance_;
    }

    SquareSubMatrix2DistMatrix() { }

    template <class T>
    void convert(const SquareSubMatrix<T>& lmat,
        dist_matrix::SparseDistMatrix<T>& dst,
        const double tol = tol_mat_elements) const;

    template <class T>
    void accumulate(const SquareSubMatrix<T>& lmat,
        dist_matrix::DistMatrix<T>& dst,
        const double tol = tol_mat_elements) const;
};

#endif
