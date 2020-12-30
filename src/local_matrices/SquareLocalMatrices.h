// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_SQUARELOCALMATRICES_H
#define MGMOL_SQUARELOCALMATRICES_H

#include "LocalMatrices.h"
#ifdef HAVE_MAGMA
class ReplicatedMatrix;
#endif

#include <vector>

template <typename DataType, typename MemorySpaceType>
class SquareLocalMatrices : public LocalMatrices<DataType, MemorySpaceType>
{
public:
    SquareLocalMatrices(const int subdiv, const int m);

    void fillUpperWithLower();
    void setDiagonal2Zero();
    void transpose();

    /*!
     *  compute trace of matrix only for elements listed in ids
     */
    double computePartialTrace(const std::vector<int>& ids, const int iloc = 0);

    /*!
     * add shift to diagonal, to shift eigenvalues
     */
    void shift(const DataType);

    /*!
     * set elements to 0 for rows/cols with gids equal to -1
     */
    void applySymmetricMask(const std::vector<std::vector<int>>& gids);

#ifdef HAVE_MAGMA
    // assign values from a ReplicatedMatrix
    void assign(const ReplicatedMatrix& src, const int iloc = 0);

    void assign(const SquareLocalMatrices<DataType, MemorySpace::Host>& src);
#endif
};

#endif
