// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Laph8.cc,v 1.4 2009/02/19 00:14:18 jeanluc Exp $
#include "Laph8.h"
const double omega = 1. / 1.5;

namespace pb
{
template <class T>
void Laph8<T>::jacobi(GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W)
{
    const double scale = -omega * invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template <class T>
void Laph8<T>::jacobi(
    GridFuncVector<T>& A, const GridFuncVector<T>& B, GridFunc<T>& W)
{
    const double scale = -omega * invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template <class T>
void Laph8<T>::jacobi(
    GridFuncVector<T>& A, const GridFuncVector<T>& B, GridFuncVector<T>& W)
{
    const double scale = -omega * invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template class Laph8<double>;
template class Laph8<float>;
} // namespace pb
