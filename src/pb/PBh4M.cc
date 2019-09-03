// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: PBh4M.cc,v 1.15 2010/01/28 22:56:31 jeanluc Exp $
#include "PBh4M.h"
#include "tools.h"

namespace pb
{

using namespace std;

// constructor
template <class T>
PBh4M<T>::PBh4M(const Grid& mygrid, DielFunc<T>& myepsilon)
    : PB<T>(mygrid, myepsilon),
      pp_(PB<T>::grid_, 1, 1, 1),
      sqrt_a_(PB<T>::epsilon_),
      inv_sqrt_a_(PB<T>::epsilon_)
{
    if (mygrid.active())
    {
        if (mygrid.ghost_pt() < minNumberGhosts())
        {
            cout << "Not enough ghosts points in PBh4M::PBh4M" << endl;
            PB<T>::grid_.mype_env().globalExit(2);
        }

        // cout<<"Constructor for PBh4M"<<endl;
        sqrt_a_.sqrt_func();
        sqrt_a_.trade_boundaries();

        inv_sqrt_a_ = sqrt_a_;
        inv_sqrt_a_.inv();
        inv_sqrt_a_.trade_boundaries();

        pp_ = 0.;
        FDoper<T>::del2_4th(sqrt_a_, pp_); // -lap
        pp_ *= -1.;
        pp_ *= inv_sqrt_a_;

        pp_.trade_boundaries();

        dot_sqrt_a_ = sqrt_a_.gdot(sqrt_a_);

        PB<T>::initialized_ = true;
    }
}

// copy constructor
// All the members will refer to the same Grid object
template <class T>
PBh4M<T>::PBh4M(const PBh4M& oper)
    : PB<T>(oper),
      pp_(oper.pp_, PB<T>::grid_),
      sqrt_a_(oper.sqrt_a_, PB<T>::grid_),
      inv_sqrt_a_(oper.inv_sqrt_a_, PB<T>::grid_)
{
    // cout<<"Copy constructor for PBh4M"<<endl;
    if (PB<T>::grid_.active())
    {
        assert(PB<T>::grid_.sizeg() > 1);
    }
    // sqrt_a_.set_grid(grid_);
    // inv_sqrt_a_.set_grid(grid_);
    // pp_.set_grid(grid_);
    dot_sqrt_a_ = oper.dot_sqrt_a_;
    // cout<<"Copy constructor for PBh4M: pp_ has "<<pp_.grid().sizeg()<<"
    // points"<<flush<<endl;
}

// construct a coarse grid operator
template <class T>
PBh4M<T> PBh4M<T>::coarseOp(const Grid& mygrid)
{
    if (!mygrid.active()) return *this;

    // cout<<"construct a coarse grid"<<endl;
    Grid coarse_G = mygrid.coarse_grid();
    DielFunc<T> ecoarse(coarse_G, PB<T>::epsilon_.epsilon_max());
    GridFunc<T> ppcoarse(coarse_G, 1, 1, 1);
    PB<T>::epsilon_.restrict3D(ecoarse);
    pp_.restrict3D(ppcoarse);

    // cout<<"construct a coarse grid operator"<<endl;
    PBh4M A(coarse_G, ecoarse, ppcoarse);
    // cout<<"New PBh4M coarse grid operator: My Grid has
    // "<<A.pp_.grid().sizeg()<<" points"<<endl; A.epsilon_.set_grid(A.grid());

    return A;
}
template <class T>
PBh4M<T> PBh4M<T>::replicatedOp(const Grid& replicated_grid)
{
    if (!PB<T>::grid_.active()) return *this;

    T* replicated_func = new T[replicated_grid.gsize()];

    this->epsilon_.init_vect(replicated_func, 'g');

    DielFunc<T> replicated_epsilon(replicated_func, replicated_grid, 0);

    this->pp_.init_vect(replicated_func, 'g');
    GridFunc<T> replicated_pp(replicated_func, replicated_grid, 1, 1, 1);

    // cout<<"construct a coarse grid operator"<<endl;
    PBh4M A(replicated_grid, replicated_epsilon, replicated_pp);
    // cout<<"New PBh4M coarse grid operator: My Grid has
    // "<<A.pp_.grid().sizeg()<<" points"<<endl; A.epsilon_.set_grid(A.grid());

    return A;
}
template <class T>
void PBh4M<T>::jacobi(GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W)
{
    if (!A.grid().active()) return;

    apply(A, W);
    W -= B; // get residual

    const double c0
        = (30. / 12.)
          * (PB<T>::inv_h2(0) + PB<T>::inv_h2(1) + PB<T>::inv_h2(2));

    const T* const vv = W.uu();
    const T* const ee = pp_.uu();
    T* aa             = A.uu();

    const int n = A.grid().sizeg();
    assert(n == static_cast<int>(B.grid().sizeg()));
    assert(n == static_cast<int>(W.grid().sizeg()));
    assert(n == static_cast<int>(pp_.grid().sizeg()));

    for (int j = 0; j < n; j++)
    {
#ifndef NDEBUG
        if (c0 + ee[j] <= 0.)
        {
            cout << "c0+ee[" << j << "]=" << c0 + ee[j] << endl;
            exit(1);
        }
#endif
        aa[j] -= (T)((0.75 / (c0 + (double)ee[j])) * (double)vv[j]);
    }

    A.set_updated_boundaries(0);
    W.set_updated_boundaries(0);
}

template <class T>
void PBh4M<T>::get_vepsilon(
    GridFunc<T>& vh, GridFunc<T>& rhod, GridFunc<T>& vepsilon)
{
    if (!vh.grid().active()) return;

#ifndef NDEBUG
    if (vh.grid().mype_env().mytask() == 0)
    {
        cout << "get_vepsilon " << endl;
    }
#endif
    PB<T>::work2_ = 0.;
    for (int i = 0; i < 3; i++)
    {
        FDoper<T>::del1_4th(vh, PB<T>::work1_, i);
        PB<T>::work2_.add_prod(PB<T>::work1_, PB<T>::work1_);
    }

    PB<T>::epsilon_.Gdepsilon_rho(rhod, vepsilon);
    vepsilon *= PB<T>::work2_;
    vepsilon *= (-1. / (8. * M_PI));
}
template <class T>
void PBh4M<T>::init(GridFunc<T>& gf_rhod)
{
    if (!gf_rhod.grid().active()) return;

    PB<T>::epsilon_.Gepsilon_rho(gf_rhod);
    sqrt_a_ = PB<T>::epsilon_;

    sqrt_a_.sqrt_func();
    sqrt_a_.trade_boundaries();

    inv_sqrt_a_ = sqrt_a_;
    inv_sqrt_a_.inv();
    inv_sqrt_a_.trade_boundaries();

    pp_ = 0.;
    FDoper<T>::del2_4th(sqrt_a_, pp_); // -lap
    pp_ *= -1.;
    pp_ *= inv_sqrt_a_;

    pp_.trade_boundaries();

    dot_sqrt_a_ = sqrt_a_.gdot(sqrt_a_);

    PB<T>::initialized_ = true;
}
template class PBh4M<double>;
template class PBh4M<float>;
} // namespace pb
