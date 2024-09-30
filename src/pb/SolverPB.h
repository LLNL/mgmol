// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: SolverPB.h,v 1.11 2010/01/28 22:56:47 jeanluc Exp $
#ifndef SOLVERPB_H
#define SOLVERPB_H

#include "Solver.h"

namespace pb
{

template <class T, typename T2>
class SolverPB : public Solver<T2>
{

private:
    T oper_;
    short nu1_;
    short nu2_;
    short max_sweeps_;
    double tol_;
    short max_nlevels_;
    bool gather_coarse_level_;

    short nb_sweeps_;
    double final_residual_;
    double final_relative_residual_;
    double residual_reduction_;

public:
    SolverPB(T& oper, const short px, const short py, const short pz)
        : Solver<T2>(px, py, pz), oper_(oper)
    {
        nu1_                 = 2; // default
        nu2_                 = 2; // default
        max_sweeps_          = 10;
        tol_                 = 1.e-16;
        max_nlevels_         = 10;
        gather_coarse_level_ = true;

        nb_sweeps_               = 0;
        final_residual_          = -1.;
        final_relative_residual_ = -1.;
        residual_reduction_      = -1.;
    };

    void setup(const short nu1, const short nu2, const short max_sweeps,
        const double tol, const short max_nlevels,
        const bool gather_coarse_level = true) override
    {
        nu1_                 = nu1;
        nu2_                 = nu2;
        max_sweeps_          = max_sweeps;
        tol_                 = tol;
        max_nlevels_         = max_nlevels;
        gather_coarse_level_ = gather_coarse_level;
    }

    bool solve(T2* phi, T2* rhs, T2* rhod, T2* vks, const char dis);
    bool solve(GridFunc<T2>& gf_phi, const GridFunc<T2>& gf_rhs,
        GridFunc<T2>& gf_rhod, GridFunc<T2>& gf_vks);
    bool solve(GridFunc<T2>& gf_phi, const GridFunc<T2>& gf_rhs) override;

    ~SolverPB() override{};

    short getNbSweeps() const override { return nb_sweeps_; }
    double getFinalResidual() const override { return final_residual_; }
    double getFinalRelativeResidual() const override
    {
        return final_relative_residual_;
    }
    double getResidualReduction() const override { return residual_reduction_; }
};

} // namespace pb

#endif
