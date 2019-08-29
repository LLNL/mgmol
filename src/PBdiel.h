// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef included_PBdiel
#define included_PBdiel

#include "MPIdata.h"
#include "Poisson.h"

#include "PBh2.h"
#include "PBh4.h"
#include "PBh4M.h"
#include "PBh4MP.h"
#include "PBh6.h"
#include "PBh8.h"
#include "SolverPB.h"

template <class T>
class PBdiel : public Poisson
{
private:
    pb::GridFunc<RHODTYPE>* rhod_;

    pb::SolverPB<T, POTDTYPE>* poisson_solver_;

public:
    // Constructor
    PBdiel(pb::Grid& grid, const short bc[3], const double e0,
        const double rho0, const double drho0)
        : Poisson(grid, bc)
    {
        T oper(Poisson::grid_, e0, rho0, drho0);
        Poisson::vepsilon_
            = new pb::GridFunc<POTDTYPE>(Poisson::grid_, bc[0], bc[1], bc[2]);
        rhod_ = nullptr;
        poisson_solver_
            = new pb::SolverPB<T, POTDTYPE>(oper, bc[0], bc[1], bc[2]);
    };

    // Destructor
    ~PBdiel() override
    {
        delete Poisson::vepsilon_;
        delete poisson_solver_;
    }

    void set_rhod(pb::GridFunc<RHODTYPE>*) override;
    void set_vh(pb::GridFunc<POTDTYPE>&);

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
