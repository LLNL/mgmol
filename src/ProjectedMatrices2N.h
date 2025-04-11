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

template <class MatrixType>
class ProjectedMatrices2N : public ProjectedMatrices<MatrixType>
{
    int bdim_;
    MatrixType* work2N_;

public:
    ProjectedMatrices2N(
        const int ndim, const bool with_spin, const double width);

    ~ProjectedMatrices2N() override;

    void assignBlocksH(MatrixType&, MatrixType&, MatrixType&, MatrixType&);

    void iterativeUpdateDMwithEigenstates(
        const double occ_width, const bool flag_reduce_T = true);
    void diagonalizeDM(std::vector<double>& occ, MatrixType& vect)
    {
        // we are assuming Gram matrix=identity
        ProjectedMatrices<MatrixType>::dm_->diagonalize('v', occ, vect);
    }
};

#endif
