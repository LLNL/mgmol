// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef REPLICATED_WORKSPACE_H
#define REPLICATED_WORKSPACE_H

#include "DistMatrix.h"
#include "Timer.h"

#include <cassert>

#define MAX_SIZE (400000)

// see S. Meyers, "More Effective C++", Item 26
// (limiting the number of objects of a class)
class ReplicatedWorkSpace
{
    DISTMATDTYPE* square_matrix_;
    int ndim_;

    DISTMATDTYPE work_double[MAX_SIZE];

    static Timer mpisum_tm_;

    ReplicatedWorkSpace()
    {
        ndim_          = 0;
        square_matrix_ = 0;
    }
    ReplicatedWorkSpace(const ReplicatedWorkSpace&);

    ~ReplicatedWorkSpace() { delete[] square_matrix_; }

public:
    Timer mpisum_tm() { return mpisum_tm_; }

    static ReplicatedWorkSpace& instance();

    void setup(const int ndim)
    {
        ndim_          = ndim;
        square_matrix_ = new DISTMATDTYPE[ndim_ * ndim_];
    }

    DISTMATDTYPE* square_matrix()
    {
        assert(square_matrix_ != 0);
        return square_matrix_;
    }

    void resetmem()
    {
        assert(square_matrix_ != 0);
        memset(square_matrix_, 0, ndim_ * ndim_ * sizeof(DISTMATDTYPE));
    }

    void mpiBcastSquareMatrix();

    void setUpperTriangularSquareMatrixToZero();

    void initSquareMatrix(const dist_matrix::DistMatrix<DISTMATDTYPE>& distmat)
    {
        distmat.matgather(square_matrix_, ndim_);
    }
    int getDim() { return ndim_; }
};

#endif
