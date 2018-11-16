// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef DFTSOLVER_H
#define DFTSOLVER_H

#include "Timer.h"
#include <iostream>

class LocGridOrbitals;
class Ions;
class MGmol;
class OrbitalsStepper;
class Energy;
class Electrostatic;
class Hamiltonian;
class ProjectedMatricesInterface;
class Rho;
class DMStrategy;

class DFTsolver
{
private:
    static Timer solve_tm_;
    static int it_scf_;

    MGmol* mgmol_strategy_;

    Hamiltonian* hamiltonian_;
    ProjectedMatricesInterface* proj_matrices_;
    Energy* energy_;
    Electrostatic* electrostat_;
    Ions& ions_;
    Rho* rho_;
    DMStrategy* dm_strategy_;

    OrbitalsStepper* orbitals_stepper_;

    // flag used to turn ON acceleration algorithm (if available) close to
    // convergence
    bool accelerate_;

    std::ostream& os_;

    double eks_history_[2];

    double sum_eig_[2];
    double deig_;
    double deig2_;
    double de_;
    double de2_;

    void printEnergy(const short) const;
    int checkConvergenceEnergy(const short step, const short max_steps);
    double evaluateEnergy(const LocGridOrbitals& orbitals, const bool flag);
    void incInnerIt() { it_scf_++; }
    bool checkPrintResidual(const short step) const;
    void dielON();
    bool testUpdatePot() const;
    bool checkConvPot() const;

public:
    DFTsolver(Hamiltonian* hamiltonian,
        ProjectedMatricesInterface* proj_matrices, Energy* energy,
        Electrostatic* electrostat, MGmol* mgmol_strategy, Ions& ions, Rho* rho,
        DMStrategy* dm_strategy, std::ostream& os);

    ~DFTsolver();

    int solve(LocGridOrbitals& orbitals, LocGridOrbitals& work_orbitals,
        Ions& ions, const short max_steps, const short iprint,
        double& last_eks);

    static void resetItCount() { it_scf_ = 0; }
    static void setItCountLarge() { it_scf_ = 1000; }

    static void printTimers(std::ostream& os);
};

#endif
