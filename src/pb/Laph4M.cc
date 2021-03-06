// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Laph4M.h"

namespace pb
{
template <class T>
void Laph4M<T>::jacobi(GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W)
{
    const double scale = -1. * jacobiFactor();

    Lap<T>::jacobi(A, B, W, scale);
}

template class Laph4M<double>;
template class Laph4M<float>;
} // namespace pb
