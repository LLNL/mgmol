// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Laph2.cc,v 1.5 2009/02/10 16:58:41 jeanluc Exp $
#include "Laph2.h"

// optimal factor according to "Multigrid" by Trottenberg, Osterlee, Schueler
// p. 73
const double omega = 6. / 7.;

namespace pb
{
template <class T>
void Laph2<T>::jacobi(GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W)
{
    const double scale = -omega * invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template <class T>
void Laph2<T>::jacobi(
    GridFuncVector<T>& A, const GridFuncVector<T>& B, GridFunc<T>& W)
{
    const double scale = -omega * invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template <class T>
void Laph2<T>::jacobi(
    GridFuncVector<T>& A, const GridFuncVector<T>& B, GridFuncVector<T>& W)
{
    const double scale = -omega * invDiagEl_;

    Lap<T>::jacobi(A, B, W, scale);
}
template class Laph2<double>;
template class Laph2<float>;
} // namespace pb
