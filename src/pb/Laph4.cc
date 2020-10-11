// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Laph4.cc,v 1.5 2008/12/10 02:08:36 jeanluc Exp $
#include "Laph4.h"
const double omega = 1. / 1.5;

namespace pb
{
template <class T>
void Laph4<T>::jacobi(GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W)
{
    const double scale = -omega * invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template <class T>
void Laph4<T>::jacobi(GridFuncVector<T, memory_space_type>& A,
    const GridFuncVector<T, memory_space_type>& B, GridFunc<T>& W)
{
    const double scale = -omega * invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template <class T>
void Laph4<T>::jacobi(GridFuncVector<T, memory_space_type>& A,
    const GridFuncVector<T, memory_space_type>& B,
    GridFuncVector<T, memory_space_type>& W)
{
    const double scale = -omega * invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template class Laph4<double>;
template class Laph4<float>;
} // namespace pb
