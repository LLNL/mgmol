// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SolverPB.h"
#include "Mgm.h"
#include "PBh2.h"
#include "PBh4.h"
#include "PBh4MP.h"
#include "PBh6.h"
#include "PBh8.h"

namespace pb
{
/*
template bool Mgm(PBh4MP&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&, const
short, const short, const double, const short, const short, const bool,
                    double&,double&,double&,short&);
template bool Mgm(PBh4M&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&, const
short, const short, const double, const short, const short, const bool,
                    double&,double&,double&,short&);
template bool Mgm(PBh2&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&, const
short, const short, const double, const short, const short, const bool,
                    double&,double&,double&,short&);
template bool Mgm(PBh4&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&, const
short, const short, const double, const short, const short, const bool,
                    double&,double&,double&,short&);
template bool Mgm(PBh6&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&, const
short, const short, const double, const short, const short, const bool,
                    double&,double&,double&,short&);
template bool Mgm(PBh8&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE>&, const
short, const short, const double, const short, const short, const bool,
                  double&,double&,double&,short&);
*/
/*
template int Vcycle(PBh4MP&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE> &,
const short, const short, const short, const bool); template int Vcycle(PBh4M&,
GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE> &, const short, const short, const
short, const bool); template int Vcycle(PBh2&, GridFunc<GFDTYPE>&, const
GridFunc<GFDTYPE> &, const short, const short, const short, const bool);
template int Vcycle(PBh4&, GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE> &, const
short, const short, const short, const bool); template int Vcycle(PBh6&,
GridFunc<GFDTYPE>&, const GridFunc<GFDTYPE> &, const short, const short, const
short, const bool); template int Vcycle(PBh8&, GridFunc<GFDTYPE>&, const
GridFunc<GFDTYPE> &, const short, const short, const short, const bool);
*/

template <class T, typename T2>
bool SolverPB<T, T2>::solve(T2* phi, T2* rhs, T2* rhod, T2* vks, const char dis)
{
    GridFunc<T2> gf_phi(oper_.grid(), Solver<T2>::bc_[0], Solver<T2>::bc_[1],
        Solver<T2>::bc_[2]);
    gf_phi.assign(phi, dis);

    GridFunc<T2> gf_work(oper_.grid(), Solver<T2>::bc_[0], Solver<T2>::bc_[1],
        Solver<T2>::bc_[2]);
    gf_work.assign(rhs, dis);

    GridFunc<T2> gf_rhod(oper_.grid(), Solver<T2>::bc_[0], Solver<T2>::bc_[1],
        Solver<T2>::bc_[2]);
    gf_rhod.assign(rhod, dis);

    oper_.init(gf_rhod);

    bool conv = Mgm(oper_, gf_phi, gf_work, max_nlevels_, max_sweeps_, tol_,
        nu1_, nu2_, gather_coarse_level_, final_residual_,
        final_relative_residual_, residual_reduction_, nb_sweeps_);

    // compute additional KS potential due to epsilon(rho)
    oper_.get_vepsilon(gf_phi, gf_rhod, gf_work);

    // convert gf_work into a double*
    gf_work.init_vect(vks, dis);

    // convert gf_phi back into a double*
    gf_phi.init_vect(phi, dis);

    return conv;
}

template <class T, typename T2>
bool SolverPB<T, T2>::solve(GridFunc<T2>& gf_phi, const GridFunc<T2>& gf_rhs,
    GridFunc<T2>& gf_rhod, GridFunc<T2>& gf_vks)
{
    oper_.init(gf_rhod);

    bool conv = Mgm(oper_, gf_phi, gf_rhs, max_nlevels_, max_sweeps_, tol_,
        nu1_, nu2_, gather_coarse_level_, final_residual_,
        final_relative_residual_, residual_reduction_, nb_sweeps_);

    // compute additional KS potential due to epsilon(rho)
    oper_.get_vepsilon(gf_phi, gf_rhod, gf_vks);

    return conv;
}

template <class T, typename T2>
bool SolverPB<T, T2>::solve(GridFunc<T2>& gf_phi, const GridFunc<T2>& gf_rhs)
{
    if (!oper_.initialized())
    {
        std::cout << "Error in SolverPB<T>::solve: operator not initialized"
                  << std::endl;
        return 0.;
    }

    bool conv = Mgm(oper_, gf_phi, gf_rhs, max_nlevels_, max_sweeps_, tol_,
        nu1_, nu2_, gather_coarse_level_, final_residual_,
        final_relative_residual_, residual_reduction_, nb_sweeps_);

    return conv;
}

// explicit instantiation declaration
#if 1
template class SolverPB<PBh4MP<double>, double>;
template class SolverPB<PBh4M<double>, double>;
template class SolverPB<PBh4<double>, double>;
template class SolverPB<PBh2<double>, double>;
template class SolverPB<PBh6<double>, double>;
template class SolverPB<PBh8<double>, double>;
template class SolverPB<PBh4MP<float>, float>;
template class SolverPB<PBh4M<float>, float>;
template class SolverPB<PBh4<float>, float>;
template class SolverPB<PBh2<float>, float>;
template class SolverPB<PBh6<float>, float>;
template class SolverPB<PBh8<float>, float>;
#else
// double
template bool SolverPB<PBh4M<double>, double>::solve(double*, double*, double*,
    const double, const double, double*, const char, const int, const int);
template bool SolverPB<PBh2<double>, double>::solve(double*, double*, double*,
    const double, const double, double*, const char, const int, const int);
template bool SolverPB<PBh4<double>, double>::solve(double*, double*, double*,
    const double, const double, double*, const char, const int, const int);
template bool SolverPB<PBh6<double>, double>::solve(double*, double*, double*,
    const double, const double, double*, const char, const int, const int);

template bool SolverPB<PBh4M<double>, float>::solve(
    GridFunc<double>&, GridFunc<double>&, const int);
template bool SolverPB<PBh2<double>, float>::solve(
    GridFunc<double>&, GridFunc<double>&, const int);
template bool SolverPB<PBh4<double>, float>::solve(
    GridFunc<double>&, GridFunc<double>&, const int);
template bool SolverPB<PBh6<double>, float>::solve(
    GridFunc<double>&, GridFunc<double>&, const int);
// float
template bool SolverPB<PBh4M<float>, float>::solve(float*, float*, float*,
    const double, const double, float*, const char, const int, const int);
template bool SolverPB<PBh2<float>, float>::solve(float*, float*, double*,
    const double, const double, float*, const char, const int, const int);
template bool SolverPB<PBh4<float>, float>::solve(float*, float*, double*,
    const double, const double, float*, const char, const int, const int);
template bool SolverPB<PBh6<float>, float>::solve(float*, float*, float*,
    const double, const double, float*, const char, const int, const int);

template bool SolverPB<PBh4M<float>, float>::solve(
    GridFunc<float>&, GridFunc<float>&, const int);
template bool SolverPB<PBh2<float>, float>::solve(
    GridFunc<float>&, GridFunc<float>&, const int);
template bool SolverPB<PBh4<float>, float>::solve(
    GridFunc<float>&, GridFunc<float>&, const int);
template bool SolverPB<PBh6<float>, float>::solve(
    GridFunc<float>&, GridFunc<float>&, const int);
#endif

} // namespace pb
