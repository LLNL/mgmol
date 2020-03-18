// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef included_ShiftedHartree
#define included_ShiftedHartree

#include "Poisson.h"

#include "ShiftedLaph4M.h"
#include "SolverLap.h"

template <class T>
class ShiftedHartree : public Poisson
{
private:
    double shift_;

    pb::SolverLap<T, POTDTYPE>* poisson_solver_;

public:
    // Constructor
    ShiftedHartree(const pb::Grid& grid, const short bc[3], const double shift)
        : Poisson(grid, bc)
    {
        shift_ = shift;
        T oper(Poisson::grid_, shift);
        poisson_solver_
            = new pb::SolverLap<T, POTDTYPE>(oper, bc[0], bc[1], bc[2]);
    };

    // Destructor
    ~ShiftedHartree() override { delete poisson_solver_; }

    void setup(const short nu1, const short nu2, const short max_sweeps,
        const double tol, const short max_nlevels,
        const bool gather_coarse_level) override
    {
        poisson_solver_->setup(
            nu1, nu2, max_sweeps, tol, max_nlevels, gather_coarse_level);
    }

    void solve(const pb::GridFunc<RHODTYPE>& rho,
        const pb::GridFunc<RHODTYPE>& rhoc) override;
};

#endif
