// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_REPLICATED_WORKSPACE_H
#define MGMOL_REPLICATED_WORKSPACE_H

#include "ReplicatedMatrix.h"
#include "Timer.h"

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrix.h"
#endif

#include <cassert>

class ReplicatedMatrix;

#define MAX_SIZE (400000)

// see S. Meyers, "More Effective C++", Item 26
// (limiting the number of objects of a class)

template <class ScalarType>
class ReplicatedWorkSpace
{
    ScalarType* square_matrix_;
    int ndim_;

    ScalarType work_double[MAX_SIZE];

    static Timer mpisum_tm_;

    ReplicatedWorkSpace()
    {
        ndim_          = 0;
        square_matrix_ = nullptr;
    }
    ReplicatedWorkSpace(const ReplicatedWorkSpace&);

    ~ReplicatedWorkSpace()
    {
        delete[] square_matrix_;
        square_matrix_ = nullptr;
    }

public:
    Timer mpisum_tm() { return mpisum_tm_; }

    static ReplicatedWorkSpace& instance();

    void setup(const int ndim)
    {
        ndim_          = ndim;
        square_matrix_ = new ScalarType[ndim_ * ndim_];
    }

    ScalarType* square_matrix()
    {
        assert(square_matrix_ != 0);
        return square_matrix_;
    }

    void resetmem()
    {
        assert(square_matrix_ != 0);
        memset(square_matrix_, 0, ndim_ * ndim_ * sizeof(ScalarType));
    }

    void mpiBcastSquareMatrix();

    void setUpperTriangularSquareMatrixToZero();

#ifdef MGMOL_USE_SCALAPACK
    void initSquareMatrix(const dist_matrix::DistMatrix<ScalarType>& tmat);
#endif
    void initSquareMatrix(const ReplicatedMatrix& mat);

    int getDim() { return ndim_; }
};

#endif
