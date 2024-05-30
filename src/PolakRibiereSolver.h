// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_PolakRibiereSolver_H
#define MGMOL_PolakRibiereSolver_H

#include "DMStrategy.h"
#include "Energy.h"
#include "Hamiltonian.h"
#include "MGmol.h"
#include "Rho.h"
#include "Timer.h"
#include <iostream>

class Ions;
class Electrostatic;
class ProjectedMatricesInterface;

template <class OrbitalsType>
class PolakRibiereSolver
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

    std::ostream& os_;

    double eks_history_[2];

    double alpha_;
    double sigma_a_;
    double sigma_b_;

    bool with_preconditioner_;

    /*!
     * residual (gradient with negative sign)
     */
    OrbitalsType* r_k_;
    OrbitalsType* r_km1_;

    /*!
     * preconditioned residual
     */
    OrbitalsType* z_k_;
    OrbitalsType* z_km1_;

    OrbitalsType* p_k_;

    double sum_eig_[2];
    double deig_;
    double deig2_;
    double de_;
    double de2_;

    void printEnergy(const short) const;
    int checkConvergenceEnergy(const short step, const short max_steps);
    double evaluateEnergy(const OrbitalsType& orbitals, const bool flag);
    void incInnerIt() { it_scf_++; }
    bool checkPrintResidual(const short step) const;
    void dielON();
    bool testUpdatePot() const;
    bool checkConvPot() const;

    bool checkWolfeConditions(
        const double trial_step_energy, const double alpha) const;

    double computeBeta(OrbitalsType& work_orbitals) const;

public:
    PolakRibiereSolver(Hamiltonian<OrbitalsType>* hamiltonian,
        ProjectedMatricesInterface* proj_matrices, Energy<OrbitalsType>* energy,
        Electrostatic* electrostat, MGmol<OrbitalsType>* mgmol_strategy,
        Ions& ions, Rho<OrbitalsType>* rho,
        DMStrategy<OrbitalsType>* dm_strategy, std::ostream& os);

    ~PolakRibiereSolver();

    int solve(OrbitalsType& orbitals, OrbitalsType& work_orbitals, Ions& ions,
        const short max_steps, const short iprint, double& last_eks);

    static void resetItCount() { it_scf_ = 0; }
    static void setItCountLarge() { it_scf_ = 1000; }

    static void printTimers(std::ostream& os);
};

#endif
