// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: SolverLap.cc,v 1.15 2010/01/28 22:56:31 jeanluc Exp $

#include "SolverLap.h"
#include "Laph2.h"
#include "Laph4.h"
#include "Laph4M.h"
#include "Laph4MP.h"
#include "Laph6.h"
#include "Laph8.h"
#include "Mgm.h"
#include "ShiftedLaph4M.h"

Timer vcycle_repwrk_tm("vcycle_repwrk");
Timer vcycle_repinit_tm("vcycle_repinit");
Timer vcycle_repend_tm("vcycle_repend");
Timer vcycle_repvcycle_tm("vcycle_repvcycle");

namespace pb
{

// explicit instantiation declaration
#if 1
// double
template class SolverLap<ShiftedLaph4M<double>, double>;
template class SolverLap<Laph4MP<double>, double>;
template class SolverLap<Laph4M<double>, double>;
template class SolverLap<Laph4<double>, double>;
template class SolverLap<Laph2<double>, double>;
template class SolverLap<Laph6<double>, double>;
template class SolverLap<Laph8<double>, double>;
// float
template class SolverLap<ShiftedLaph4M<float>, float>;
template class SolverLap<Laph4MP<float>, float>;
template class SolverLap<Laph4M<float>, float>;
template class SolverLap<Laph4<float>, float>;
template class SolverLap<Laph2<float>, float>;
template class SolverLap<Laph6<float>, float>;
template class SolverLap<Laph8<float>, float>;
#else
// double
template bool SolverLap<Laph4M<double>, double>::solve(
    double*, double*, const char, const int, const int);
template bool SolverLap<Laph2<double>, double>::solve(
    double*, double*, const char, const int, const int);
template bool SolverLap<Laph4<double>, double>::solve(
    double*, double*, const char, const int, const int);
template bool SolverLap<Laph6<double>, double>::solve(
    double*, double*, const char, const int, const int);

template bool SolverLap<Laph4M<double>, double>::solve(
    GridFunc<double>&, GridFunc<double>&, const int);
template bool SolverLap<Laph2<double>, double>::solve(
    GridFunc<double>&, GridFunc<double>&, const int);
template bool SolverLap<Laph4<double>, double>::solve(
    GridFunc<double>&, GridFunc<double>&, const int);
template bool SolverLap<Laph6<double>, double>::solve(
    GridFunc<double>&, GridFunc<double>&, const int);
// float
template bool SolverLap<Laph4M<float>, float>::solve(
    float*, float*, const char, const int, const int);
template bool SolverLap<Laph2<float>, float>::solve(
    float*, float*, const char, const int, const int);
template bool SolverLap<Laph4<float>, float>::solve(
    float*, float*, const char, const int, const int);
template bool SolverLap<Laph6<float>, float>::solve(
    float*, float*, const char, const int, const int);

template bool SolverLap<Laph4M<float>, float>::solve(
    GridFunc<float>&, GridFunc<float>&, const int);
template bool SolverLap<Laph2<float>, float>::solve(
    GridFunc<float>&, GridFunc<float>&, const int);
template bool SolverLap<Laph4<float>, float>::solve(
    GridFunc<float>&, GridFunc<float>&, const int);
template bool SolverLap<Laph6<float>, float>::solve(
    GridFunc<float>&, GridFunc<float>&, const int);
#endif
/*
template bool Mgm(Laph4MP&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&, const
short, const short, const double, const short, const short, const bool,
                    double&,double&,double&,short&);
template bool Mgm(Laph4M&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&, const
short, const short, const double, const short, const short, const bool,
                    double&,double&,double&,short&);
template bool Mgm(Laph2&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&, const
short, const short, const double, const short, const short, const bool,
                    double&,double&,double&,short&);
template bool Mgm(Laph4&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&, const
short, const short, const double, const short, const short, const bool,
                    double&,double&,double&,short&);
template bool Mgm(Laph6&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&, const
short, const short, const double, const short, const short, const bool,
                    double&,double&,double&,short&);
template bool Mgm(ShiftedLaph4M&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&,
const short, const short, const double, const short, const short, const bool,
                    double&,double&,double&,short&);
*/
/*
template int Vcycle(ShiftedLaph4M&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>
&, const short, const short, const short, const bool); template int
Vcycle(Laph4MP&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE> &, const short,
const short, const short, const bool); template int Vcycle(Laph4M&,
GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE> &, const short, const short, const
short, const bool); template int Vcycle(Laph2&, GridFunc<GFDTYPE>&, const
GridFunc<GFDTYPE> &, const short, const short, const short, const bool);
template int Vcycle(Laph4&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE> &, const
short, const short, const short, const bool); template int Vcycle(Laph6&,
GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE> &, const short, const short, const
short, const bool);
*/
template <class T, typename T2>
bool SolverLap<T, T2>::solve(T2* phi, T2* rhs, const char dis)
{
    GridFunc<T2> gf_phi(oper_.grid(), Solver<T2>::bc_[0], Solver<T2>::bc_[1],
        Solver<T2>::bc_[2]);
    gf_phi.assign(phi, dis);
    GridFunc<T2> gf_work(oper_.grid(), Solver<T2>::bc_[0], Solver<T2>::bc_[1],
        Solver<T2>::bc_[2]);
    gf_work.assign(rhs, dis);

    bool conv = Mgm(oper_, gf_phi, gf_work, max_nlevels_, max_sweeps_, tol_,
        nu1_, nu2_, gather_coarse_level_, final_residual_,
        final_relative_residual_, residual_reduction_, nb_sweeps_);

    if (Solver<T2>::fully_periodic_) gf_phi.average0();

    // convert gf_phi back into a double*
    gf_phi.init_vect(phi, dis);

    return conv;
}

template <class T, typename T2>
bool SolverLap<T, T2>::solve(GridFunc<T2>& gf_phi, GridFunc<T2>& gf_rhs)
{
    bool conv = Mgm(oper_, gf_phi, gf_rhs, max_nlevels_, max_sweeps_, tol_,
        nu1_, nu2_, gather_coarse_level_, final_residual_,
        final_relative_residual_, residual_reduction_, nb_sweeps_);

    if (Solver<T2>::fully_periodic_) gf_phi.average0();

    return conv;
}

} // namespace pb
