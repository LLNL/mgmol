// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ELECTROSTATIC_H
#define MGMOL_ELECTROSTATIC_H

#include "GridFunc.h"
#include "Poisson.h"
#include "Rho.h"
#include "Timer.h"

class Ions;
class LocGridOrbitals;
class Potentials;

class Electrostatic
{
    short laptype_;
    short bc_[3];
    pb::Grid* pbGrid_;
    bool diel_flag_;

    Poisson* poisson_solver_;

    double Evh_rho_;
    double Evh_rhoc_;
    double Evhold_rho_;
    double eepsilon_;

    pb::GridFunc<RHODTYPE>* grhoc_;
    pb::GridFunc<RHODTYPE>* grhod_;

    int iterative_index_;

    static Timer solve_tm_;

public:
    Electrostatic(const short lap_type, const short bc[3],
        const double screening_const = 0.);
    ~Electrostatic();
    static Timer solve_tm() { return solve_tm_; }

    void setup(const short max_sweeps);
    void setupPB(const double rho0, const double drho0, Potentials& pot);

    void setupRhoc(RHODTYPE* rhoc);
    void fillFuncAroundIons(const Ions& ions);
    void computeVh(const Ions& ions, Rho<LocGridOrbitals>& rho, Potentials& pot);
    void computeVh(const pb::GridFunc<POTDTYPE>& vhinit, const Ions& ions,
        Rho<LocGridOrbitals>& rho, Potentials& pot);
    void setupInitialVh(const POTDTYPE* const);
    void setupInitialVh(const pb::GridFunc<POTDTYPE>&);
    void computeVhRho(Rho<LocGridOrbitals>& rho);
    void resetSolution() { poisson_solver_->resetVh(); }

    const pb::GridFunc<POTDTYPE>& getVh() const;

    double evhRho() const { return Evh_rho_; }
    double evhRhoc() const { return Evh_rhoc_; }
    double eEpsilon() const { return eepsilon_; }
};

#endif
