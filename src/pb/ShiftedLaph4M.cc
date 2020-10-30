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
void ShiftedLaph4M<T>::apply(GridFuncVector<T, memory_space_type>& A,
    GridFuncVector<T, memory_space_type>& B)
{
    assert(A.size() == B.size());
    const int nfunc = (int)A.size();
    A.trade_boundaries();
    GridFunc<T> C(A.getGridFunc(0).grid(), A.getGridFunc(0).bc(0),
        A.getGridFunc(0).bc(1), A.getGridFunc(0).bc(2));
    for (int k = 0; k < nfunc; k++)
    {
        FDoper<T>::del2_4th_Mehr(A.getGridFunc(k), B.getGridFunc(k));
        FDoper<T>::rhs_4th_Mehr1(A.getGridFunc(k), C);
        B.getGridFunc(k).axpy(lambda2_, C);
    }
}
template <class T>
void ShiftedLaph4M<T>::jacobi(
    GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W)
{
    const double scale = -invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template <class T>
void ShiftedLaph4M<T>::jacobi(GridFuncVector<T, memory_space_type>& A,
    const GridFuncVector<T, memory_space_type>& B, GridFunc<T>& W)
{
    const double scale = -invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template <class T>
void ShiftedLaph4M<T>::jacobi(GridFuncVector<T, memory_space_type>& A,
    const GridFuncVector<T, memory_space_type>& B,
    GridFuncVector<T, memory_space_type>& W)
{
    const double scale = -invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template class ShiftedLaph4M<double>;
template class ShiftedLaph4M<float>;
} // namespace pb
