// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef _HARTREE_CG_H_
#define _HARTREE_CG_H_

#include "PCGSolver.h"
#include "Poisson.h"

template <class T>
class Hartree_CG : public Poisson
{
private:
    PCGSolver<T, POTDTYPE>* poisson_solver_;

public:
    // Constructor
    Hartree_CG(const pb::Grid& grid, const short bc[3]) : Poisson(grid, bc)
    {
        T oper(Poisson::grid_);
        poisson_solver_ = new PCGSolver<T, POTDTYPE>(oper, bc[0], bc[1], bc[2]);
    };

    // Destructor
    ~Hartree_CG() override { delete poisson_solver_; }

    void setup(const short nu1, const short nu2, const short max_sweeps,
        const double tol, const short max_nlevels,
        const bool gather_coarse_level = true) override
    {
        (void)gather_coarse_level;
        poisson_solver_->setup(nu1, nu2, max_sweeps, tol, max_nlevels);
    }

    void solve(
        const pb::GridFunc<RHODTYPE>& rho, const pb::GridFunc<RHODTYPE>& rhoc) override;
};

#endif
