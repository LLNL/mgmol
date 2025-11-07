// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef MGMOL_MVPSOLVER_H
#define MGMOL_MVPSOLVER_H

#include "Energy.h"
#include "Hamiltonian.h"
#include "MGmol.h"
#include "Rho.h"
#include "Timer.h"

class Ions;
class Electrostatic;
template <class OrbitalsType>
class ProjectedMatrices;

template <class OrbitalsType, class MatrixType>
class MVPSolver
{
private:
    const MPI_Comm comm_;
    std::ostream& os_;

    const short n_inner_steps_;

    const bool use_old_dm_;
    Ions& ions_;

    int numst_;
    double mixing_;

    /*!
     * tolerance on energy slope in inner iterations
     */
    double tol_de0_;

    Rho<OrbitalsType>* rho_;
    Energy<OrbitalsType>* energy_;
    Electrostatic* electrostat_;
    Hamiltonian<OrbitalsType>* hamiltonian_;

    MGmol<OrbitalsType>* mgmol_strategy_;

    MatrixType* work_;
    ProjectedMatrices<MatrixType>* proj_mat_work_;

    static Timer solve_tm_;
    static Timer target_tm_;

    double evaluateDerivative(
        MatrixType& dm2Ninit, MatrixType& delta_dm, const double ts0);
    void buildTarget_MVP(MatrixType& h11, MatrixType& s11, MatrixType& target);

public:
    MVPSolver(MPI_Comm comm, std::ostream& os, Ions& ions, Rho<OrbitalsType>*,
        Energy<OrbitalsType>*, Electrostatic*, Hamiltonian<OrbitalsType>*,
        MGmol<OrbitalsType>* mgmol_strategy, const int numst, const double kbT,
        const std::vector<std::vector<int>>& global_indexes,
        const short n_inner_steps, const double mixing, const double tol_de0,
        const bool use_old_dm);
    ~MVPSolver();

    int solve(OrbitalsType& orbitals);
    static void printTimers(std::ostream& os);
};

#endif
