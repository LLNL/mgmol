// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: ShiftedLaph4M.cc,v 1.5 2009/02/19 00:14:26 jeanluc Exp $
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
void ShiftedLaph4M<T>::apply(GridFuncVector<T>& A, GridFuncVector<T>& B)
{
    assert(A.size() == B.size());
    const int nfunc = (int)A.size();
    A.trade_boundaries();
    GridFunc<T> C(
        A.func(0).grid(), A.func(0).bc(0), A.func(0).bc(1), A.func(0).bc(2));
    for (int k = 0; k < nfunc; k++)
    {
        FDoper<T>::del2_4th_Mehr(A.func(k), B.func(k));
        FDoper<T>::rhs_4th_Mehr1(A.func(k), C);
        B.func(k).axpy(lambda2_, C);
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
void ShiftedLaph4M<T>::jacobi(
    GridFuncVector<T>& A, const GridFuncVector<T>& B, GridFunc<T>& W)
{
    const double scale = -invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template <class T>
void ShiftedLaph4M<T>::jacobi(
    GridFuncVector<T>& A, const GridFuncVector<T>& B, GridFuncVector<T>& W)
{
    const double scale = -invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template class ShiftedLaph4M<double>;
template class ShiftedLaph4M<float>;
} // namespace pb
