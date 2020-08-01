// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "PBh4.h"

namespace pb
{
template <class T>
void PBh4<T>::pb_4th(GridFunc<T>& A, GridFunc<T>& B)
{
    B = 0.;

    for (int i = 0; i < 3; i++)
    {

        FDoper<T>::del1_4th(PB<T>::epsilon_, PB<T>::work2_, i);
        FDoper<T>::del1_4th(A, PB<T>::work1_, i);
        B.substract_prod(PB<T>::work1_, PB<T>::work2_);
    }

    FDoper<T>::del2_4th(A, PB<T>::work1_);
    B.add_prod(PB<T>::work1_, PB<T>::epsilon_);
}
template <class T>
void PBh4<T>::jacobi(GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W)
{
    apply(A, W);
    W -= B; // get residual

    const double c0
        = 1.5 * (30. / 12.)
          * (PB<T>::inv_h2(0) + PB<T>::inv_h2(1) + PB<T>::inv_h2(2));

    if (!PB<T>::epsilon_.updated_boundaries())
        PB<T>::epsilon_.trade_boundaries();
    assert(PB<T>::epsilon_.fully_periodic());

    A.jacobi(W, PB<T>::epsilon_, c0);

    A.set_updated_boundaries(0);
    W.set_updated_boundaries(0);
}
template <class T>
void PBh4<T>::get_vepsilon(
    GridFunc<T>& vh, GridFunc<T>& rho, GridFunc<T>& vepsilon)
{
    PB<T>::work2_ = 0.;
    for (int i = 0; i < 3; i++)
    {
        FDoper<T>::del1_4th(vh, PB<T>::work1_, i);
        PB<T>::work2_.add_prod(PB<T>::work1_, PB<T>::work1_);
    }

    PB<T>::epsilon_.Gdepsilon_rho(rho, vepsilon);
    vepsilon *= PB<T>::work2_;
    vepsilon *= (-1. / (8. * M_PI));
}
template <class T>
PBh4<T> PBh4<T>::replicatedOp(const Grid& replicated_grid)
{
    T* replicated_func = new T[replicated_grid.gsize()];

    this->epsilon_.init_vect(replicated_func, 'g');

    DielFunc<T> replicated_epsilon(replicated_grid);
    replicated_epsilon.assign(replicated_func, 0);

    delete[] replicated_func;

    PBh4 A(replicated_grid, replicated_epsilon);

    return A;
}
template class PBh4<double>;
template class PBh4<float>;
} // namespace pb
