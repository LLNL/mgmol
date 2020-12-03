// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef MGMOL_SQUARESUBMATRIX
#define MGMOL_SQUARESUBMATRIX

#include <map>
#include <vector>

template <class T>
class SquareSubMatrix
{
private:
    // number of rows and columns
    const int n_;

    // matrix elements
    std::vector<T> data_;

    // global indexes
    const std::vector<int> gids_;

    // mapping from global index to local storage
    std::map<int, int> gid2loc_;

public:
    SquareSubMatrix(const std::vector<int>& gid);

    T* data() { return data_.data(); }

    int ld() const { return n_; }

    const std::vector<int>& getGids() const { return gids_; }

    void setLocalValue(const int i, const int j, const double val)
    {
        data_[j * n_ + i] = val;
    }

    void addValue(const int gi, const int gj, const double val)
    {
        data_[gid2loc_[gj] * n_ + gid2loc_[gi]] += val;
    }

    void addLocalValue(const int i, const int j, const double val)
    {
        data_[j * n_ + i] += val;
    }

    T getLocalValue(const int i, const int j) const
    {
        return data_[j * n_ + i];
    }
};

#endif
