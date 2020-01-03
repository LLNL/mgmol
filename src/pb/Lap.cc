// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Lap.h"

namespace pb
{

template <class T>
double Lap<T>::energyES(GridFunc<T>& v, GridFunc<T>& rho)
{

    double g = 0.5 * dot(v, rho);

    if (v.mype_env().mytask() == 0)
        std::cout << " ES Energy = " << g << std::endl;

    return g;
}
template <class T>
void Lap<T>::jacobi(
    GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W, const double scale)
{
    apply(A, W);
    W -= B;

    A.axpy(scale, W);

    A.set_updated_boundaries(0);
    W.set_updated_boundaries(0);
}
template <class T>
void Lap<T>::jacobi(GridFuncVector<T>& A, const GridFuncVector<T>& B,
    GridFunc<T>& W, const double scale)
{
    assert(A.size() == B.size());

    A.trade_boundaries();

    const int nfunc = (int)A.size();
#pragma omp parallel for
    for (int k = 0; k < nfunc; k++)
    {
        apply(A.func(k), W);
        W -= (B.func(k));
        A.func(k).axpy(scale, W);
        A.func(k).set_updated_boundaries(false);
    }

    A.set_updated_boundaries(false);
}
template <class T>
void Lap<T>::jacobi(GridFuncVector<T>& A, const GridFuncVector<T>& B,
    GridFuncVector<T>& W, const double scale)
{
    assert(A.size() == B.size());

    A.trade_boundaries();

    const int nfunc = (int)A.size();
#pragma omp parallel for
    for (int k = 0; k < nfunc; k++)
    {
        apply(A.func(k), W.func(k));
        W.func(k) -= (B.func(k));
        A.func(k).axpy(scale, W.func(k));
        A.func(k).set_updated_boundaries(false);
    }

    A.set_updated_boundaries(false);
}
template class Lap<double>;
template class Lap<float>;
} // namespace pb
