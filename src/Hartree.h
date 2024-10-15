// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef included_Hartree
#define included_Hartree

#include "Poisson.h"
#include "PoissonInterface.h"

#include "SolverLap.h"

template <class T>
class Hartree : public Poisson
{
private:
    pb::SolverLap<T, POTDTYPE>* poisson_solver_;

public:
    // Constructor
    Hartree(const pb::Grid& grid, const short bc[3]) : Poisson(grid, bc)
    {
        T oper(Poisson::grid_);
        poisson_solver_
            = new pb::SolverLap<T, POTDTYPE>(oper, bc[0], bc[1], bc[2]);
    };

    // Destructor
    ~Hartree() override { delete poisson_solver_; }

    void setup(const short nu1, const short nu2, const short max_sweeps,
        const double tol, const short max_nlevels,
        const bool gather_coarse_level = true) override
    {
        poisson_solver_->setup(
            nu1, nu2, max_sweeps, tol, max_nlevels, gather_coarse_level);
    }

    void solve(const pb::GridFunc<RHODTYPE>& rho,
        const pb::GridFunc<RHODTYPE>& rhoc) override;
};

#endif
