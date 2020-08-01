// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "PBh8.h"

namespace pb
{
template <class T>
void PBh8<T>::pb_8th(GridFunc<T>& A, GridFunc<T>& B)
{
    B = 0.;
    // work1_=0.;
    // work2_=0.;

    for (int i = 0; i < 3; i++)
    {

        FDoper<T>::del1_8th(PB<T>::epsilon_, PB<T>::work2_, i);
        FDoper<T>::del1_8th(A, PB<T>::work1_, i);
        B.substract_prod(PB<T>::work1_, PB<T>::work2_);
    }

    FDoper<T>::del2_8th(A, PB<T>::work1_);
    B.add_prod(PB<T>::work1_, PB<T>::epsilon_);
}
template <class T>
void PBh8<T>::jacobi(GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W)
{
    apply(A, W);
    W -= B; // get residual

    T* u = A.uu(0);
    T* v = W.uu(0);
    T* e = PB<T>::epsilon_.uu(0);
    const double c0
        = 1.5 * (1435. / 504.)
          * (PB<T>::inv_h2(0) + PB<T>::inv_h2(1) + PB<T>::inv_h2(2));

    for (unsigned int j = 0; j < A.grid().sizeg(); j++)
    {
        assert(e[j] > 0.);
        double uval = (double)u[j] - (1. / (c0 * (double)e[j])) * (double)v[j];
        u[j]        = (T)uval;
        //	u[j] -= (T)(1./(c0*(double)e[j]))*(double)v[j];
    }

    A.set_updated_boundaries(0);
    W.set_updated_boundaries(0);
}
template <class T>
void PBh8<T>::get_vepsilon(
    GridFunc<T>& vh, GridFunc<T>& rho, GridFunc<T>& vepsilon)
{
    PB<T>::work2_ = 0.;
    for (int i = 0; i < 3; i++)
    {
        FDoper<T>::del1_8th(vh, PB<T>::work1_, i);
        PB<T>::work2_.add_prod(PB<T>::work1_, PB<T>::work1_);
    }

    PB<T>::epsilon_.Gdepsilon_rho(rho, vepsilon);
    vepsilon *= PB<T>::work2_;
    vepsilon *= (-1. / (8. * M_PI));
}
template <class T>
PBh8<T> PBh8<T>::replicatedOp(const Grid& replicated_grid)
{
    T* replicated_func = new T[replicated_grid.gsize()];

    this->epsilon_.init_vect(replicated_func, 'g');

    DielFunc<T> replicated_epsilon(replicated_grid);
    replicated_epsilon.assign(replicated_func, 0);

    delete[] replicated_func;

    PBh8 A(replicated_grid, replicated_epsilon);

    return A;
}
template class PBh8<double>;
template class PBh8<float>;
} // namespace pb
