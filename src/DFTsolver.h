// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DFTSOLVER_H
#define MGMOL_DFTSOLVER_H

#include "DMStrategy.h"
#include "DielectricControl.h"
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

template <class OrbitalsType>
class DFTsolver
{
private:
    static Timer solve_tm_;
    static int it_scf_;

    MGmol<OrbitalsType>* mgmol_strategy_;

    Hamiltonian<OrbitalsType>* hamiltonian_;
    ProjectedMatricesInterface* proj_matrices_;
    Energy<OrbitalsType>* energy_;
    Electrostatic* electrostat_;
    Ions& ions_;
    Rho<OrbitalsType>* rho_;
    DMStrategy<OrbitalsType>* dm_strategy_;

    OrbitalsStepper<OrbitalsType>* orbitals_stepper_;

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

    DielectricControl diel_control_;

    void printEnergy(const short) const;
    int checkConvergenceEnergy(const short step, const short max_steps);
    double evaluateEnergy(const OrbitalsType& orbitals, const bool flag);
    void incInnerIt() { it_scf_++; }
    bool checkPrintResidual(const short step) const;
    bool testUpdatePot() const;
    bool checkConvPot() const;

public:
    DFTsolver(Hamiltonian<OrbitalsType>* hamiltonian,
        ProjectedMatricesInterface* proj_matrices, Energy<OrbitalsType>* energy,
        Electrostatic* electrostat, MGmol<OrbitalsType>* mgmol_strategy,
        Ions& ions, Rho<OrbitalsType>* rho,
        DMStrategy<OrbitalsType>* dm_strategy, std::ostream& os);

    ~DFTsolver();

    int solve(OrbitalsType& orbitals, OrbitalsType& work_orbitals, Ions& ions,
        const short max_steps, const short iprint, double& last_eks);

    static void resetItCount() { it_scf_ = 0; }
    static void setItCountLarge() { it_scf_ = 1000; }

    static void printTimers(std::ostream& os);
};
// Instantiate static variables here to avoid clang warnings
template <class OrbitalsType>
Timer DFTsolver<OrbitalsType>::solve_tm_("solve");
template <class OrbitalsType>
int DFTsolver<OrbitalsType>::it_scf_ = 0;
#endif
