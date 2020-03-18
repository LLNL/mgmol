// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "SquareSubMatrix.h"

#include <string.h>

template <class T>
SquareSubMatrix<T>::SquareSubMatrix(const std::vector<int>& gids)
    : n_(gids.size()), gids_(gids)
{
    size_t nels = n_ * n_;

    data_.resize(nels);

    memset(data_.data(), 0., nels * sizeof(T));

    int i = 0;
    for (auto v : gids_)
    {
        gid2loc_.insert(std::pair<int, int>(v, i));
        i++;
    }
}

template class SquareSubMatrix<double>;
