// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "ShiftedLaph4M.h"

namespace pb
{

// A->B
template <class T>
void ShiftedLaph4M<T>::apply(GridFunc<T>& A, GridFunc<T>& B)
{
    FDoper<T>::del2_4th_Mehr(A, B);
    GridFunc<T> C(A.grid(), A.bc(0), A.bc(1), A.bc(2));
    FDoper<T>::rhs_4th_Mehr1(A, C);
    B.axpy(lambda2_, C);
    B.set_bc(A.bc(0), A.bc(1), A.bc(2));
}

template <class T>
void ShiftedLaph4M<T>::jacobi(
    GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W)
{
    const double scale = -1. * jacobiFactor();

    Lap<T>::jacobi(A, B, W, scale);
}

template class ShiftedLaph4M<double>;
template class ShiftedLaph4M<float>;
} // namespace pb
