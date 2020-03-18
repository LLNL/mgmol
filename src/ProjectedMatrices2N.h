// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_PROJECTED_MATRICES2N_H
#define MGMOL_PROJECTED_MATRICES2N_H

#include "ProjectedMatrices.h"

class ProjectedMatrices2N : public ProjectedMatrices
{
    int bdim_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* work2N_;

public:
    ProjectedMatrices2N(const int ndim, const bool with_spin);

    ~ProjectedMatrices2N() override;

    void assignBlocksH(dist_matrix::DistMatrix<DISTMATDTYPE>&,
        dist_matrix::DistMatrix<DISTMATDTYPE>&,
        dist_matrix::DistMatrix<DISTMATDTYPE>&,
        dist_matrix::DistMatrix<DISTMATDTYPE>&);
    double mu() const { return mu_; }

    void iterativeUpdateDMwithEigenstates(const double occ_width, const int nel,
        const int iterative_index, const bool flag_reduce_T = true);
    void diagonalizeDM(std::vector<DISTMATDTYPE>& occ,
        dist_matrix::DistMatrix<DISTMATDTYPE>& vect)
    {
        // we are assuming Gram matrix=identity
        dm_->diagonalize('v', occ, vect);
    }
};

#endif
