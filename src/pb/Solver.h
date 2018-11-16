// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Solver.h,v 1.14 2010/01/28 22:56:47 jeanluc Exp $
#ifndef SOLVE_H
#define SOLVE_H

#include "GridFunc.h"

namespace pb
{
template <typename T>
class Solver
{

protected:
    short bc_[3];
    bool fully_periodic_;

public:
    Solver(const short px, const short py, const short pz)
    {
        bc_[0] = px;
        bc_[1] = py;
        bc_[2] = pz;

        fully_periodic_ = ((bc_[0] == 1) && (bc_[1] == 1) && (bc_[2] == 1));
    }

    virtual bool solve(GridFunc<T>&, GridFunc<T>&) = 0;

    virtual void setup(const short nu1, const short nu2, const short max_sweeps,
        const double tol, const short max_nlevels,
        const bool gather_coarse_level = true)
        = 0;

    virtual ~Solver() {}

    virtual short getNbSweeps() const               = 0;
    virtual double getFinalResidual() const         = 0;
    virtual double getFinalRelativeResidual() const = 0;
    virtual double getResidualReduction() const     = 0;
};

} // namespace pb

#endif
