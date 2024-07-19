// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef _PCG_SOLVER_H_
#define _PCG_SOLVER_H_

#include "Control.h"
#include "Lap.h"
#include "LapFactory.h"
#include "global.h"

#include <vector>

template <class T, typename T2>
class PCGSolver
{

private:
    std::vector<pb::Grid*> grid_;
    short lap_type_;
    short bc_[3];
    bool fully_periodic_;
    // operators
    T oper_;
    std::vector<pb::Lap<POISSONPRECONDTYPE>*> pc_oper_;
    std::vector<pb::GridFunc<POISSONPRECONDTYPE>*> gf_work_;
    std::vector<pb::GridFunc<POISSONPRECONDTYPE>*> gf_rcoarse_;
    std::vector<pb::GridFunc<POISSONPRECONDTYPE>*> gf_newv_;
    // solver params
    int maxiters_;
    double tol_;
    double final_residual_;
    double residual_reduction_;
    // precon params
    short nu1_;
    short nu2_;
    short max_nlevels_;
    short nlevels_;
    bool is_pc_setup_;

public:
    PCGSolver(T& oper, const short px, const short py, const short pz)
        : oper_(oper)
    {
        maxiters_           = 10; // default
        nu1_                = 2; // default
        nu2_                = 2; // default
        tol_                = 1.e-16;
        max_nlevels_        = 10;
        final_residual_     = -1.;
        residual_reduction_ = -1.;

        // boundary conditions
        bc_[0]          = px;
        bc_[1]          = py;
        bc_[2]          = pz;
        fully_periodic_ = ((bc_[0] == 1) && (bc_[1] == 1) && (bc_[2] == 1));
        // fine grid info
        //        pb::Grid* mygrid=new pb::Grid(oper.grid());
        //        grid_.push_back(mygrid);
        // fine grid operator
        Control& ct  = *(Control::instance());
        lap_type_    = ct.lap_type;
        is_pc_setup_ = false;
        //        pb::Lap* myoper = LapFactory::createLap(*grid_[0],lap_type_);
        //        pc_oper_.push_back(myoper);
    };

    void setup(const short nu1, const short nu2, const short max_sweeps,
        const double tol, const short max_nlevels)
    {
        maxiters_    = max_sweeps;
        nu1_         = nu1;
        nu2_         = nu2;
        tol_         = tol;
        max_nlevels_ = max_nlevels;
        setupPrecon();
    }

    void clear();

    void setupPrecon();

    void preconSolve(pb::GridFunc<POISSONPRECONDTYPE>& gf_v,
        const pb::GridFunc<POISSONPRECONDTYPE>& gf_f, const short level = 0);

    bool solve(pb::GridFunc<T2>& gf_phi, pb::GridFunc<T2>& gf_rhs);

    bool solve(T2* phi, T2* rhs, const char dis);

    double getFinalResidual() const { return final_residual_; }
    double getResidualReduction() const { return residual_reduction_; }

    T* getOperator() { return &oper_; }

    // Destructor
    ~PCGSolver() { clear(); }
};

#endif
