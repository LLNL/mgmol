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

#include "DistMatrix.h"
#include "Timer.h"

#include <cassert>

#define MAX_SIZE (400000)

// see S. Meyers, "More Effective C++", Item 26
// (limiting the number of objects of a class)

template <class T>
class ReplicatedWorkSpace
{
    T* square_matrix_;
    int ndim_;

    T work_double[MAX_SIZE];

    static Timer mpisum_tm_;

    ReplicatedWorkSpace()
    {
        ndim_          = 0;
        square_matrix_ = nullptr;
    }
    ReplicatedWorkSpace(const ReplicatedWorkSpace&);

    ~ReplicatedWorkSpace() { delete[] square_matrix_; }

public:
    Timer mpisum_tm() { return mpisum_tm_; }

    static ReplicatedWorkSpace& instance();

    void setup(const int ndim)
    {
        ndim_          = ndim;
        square_matrix_ = new T[ndim_ * ndim_];
    }

    T* square_matrix()
    {
        assert(square_matrix_ != 0);
        return square_matrix_;
    }

    void resetmem()
    {
        assert(square_matrix_ != 0);
        memset(square_matrix_, 0, ndim_ * ndim_ * sizeof(T));
    }

    void mpiBcastSquareMatrix();

    void setUpperTriangularSquareMatrixToZero();

    void initSquareMatrix(const dist_matrix::DistMatrix<T>& distmat)
    {
        distmat.allgather(square_matrix_, ndim_);
    }
    int getDim() { return ndim_; }
};

#endif
