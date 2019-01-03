// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DFTSOLVER_H
#define MGMOL_DFTSOLVER_H

#include "Energy.h"
#include "Hamiltonian.h"
#include "MGmol.h"
#include "OrbitalsStepper.h"
#include "Rho.h"
#include "Timer.h"
#include <iostream>

class Ions;
class Electrostatic;
class ProjectedMatricesInterface;
class DMStrategy;

template <class T>
class DFTsolver
{
private:
    static Timer solve_tm_;
    static int it_scf_;

    MGmol<T>* mgmol_strategy_;

    Hamiltonian<T>* hamiltonian_;
    ProjectedMatricesInterface* proj_matrices_;
    Energy<T>* energy_;
    Electrostatic* electrostat_;
    Ions& ions_;
    Rho<T>* rho_;
    DMStrategy* dm_strategy_;

    OrbitalsStepper<T>* orbitals_stepper_;

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
    double evaluateEnergy(const T& orbitals, const bool flag);
    void incInnerIt() { it_scf_++; }
    bool checkPrintResidual(const short step) const;
    void dielON();
    bool testUpdatePot() const;
    bool checkConvPot() const;

public:
    DFTsolver(Hamiltonian<T>* hamiltonian,
        ProjectedMatricesInterface* proj_matrices, Energy<T>* energy,
        Electrostatic* electrostat, MGmol<T>* mgmol_strategy, Ions& ions,
        Rho<T>* rho,
        DMStrategy* dm_strategy, std::ostream& os);

    ~DFTsolver();

    int solve(T& orbitals, T& work_orbitals,
        Ions& ions, const short max_steps, const short iprint,
        double& last_eks);

    static void resetItCount() { it_scf_ = 0; }
    static void setItCountLarge() { it_scf_ = 1000; }

    static void printTimers(std::ostream& os);
};

#endif
