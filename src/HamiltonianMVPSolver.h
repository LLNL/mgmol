// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_HAMILTONIANMVP_SOLVER_H_
#define MGMOL_HAMILTONIANMVP_SOLVER_H_

#include "Energy.h"
#include "MGmol.h"
#include "Rho.h"
#include "Timer.h"

class Ions;
template <class T>
class MGmol;
class Electrostatic;
template <class MatrixType>
class ProjectedMatrices2N;
template <class MatrixType>
class ProjectedMatrices;

template <class T1, class T2, class T3>
class HamiltonianMVPSolver
{

private:
    std::ostream& os_;

    short n_inner_steps_;

    Ions& ions_;

    Rho<T3>* rho_;
    Energy<T3>* energy_;
    Electrostatic* electrostat_;
    MGmol<T3>* mgmol_strategy_;

    int numst_;

    /*!
     * If this flag is on, try shortening interval for line minimization
     * until successful (interpolation coefficient larger than 0)
     */
    bool try_shorter_intervals_;

    /*!
     * "variable" matrix defining "variable" DM through diagonalization
     * keep values from previous call to solve() function
     */
    T1* hmatrix_;

    /*!
     * Initial natrix in last solve.
     * Used to reset hmatrix_ when move not accepted in outer solver
     */
    T1* initial_hmatrix_;

    static Timer solve_tm_;
    static Timer target_tm_;

public:
    HamiltonianMVPSolver(std::ostream& os, Ions& ions, Rho<T3>* rho,
        Energy<T3>* energy, Electrostatic* electrostat,
        MGmol<T3>* mgmol_strategy, const int numst, const short n_inner_steps,
        const T1& hinit, const bool try_shorter_intervals = false);
    ~HamiltonianMVPSolver();
    int solve(T3& orbitals);
    void reset();
    void printTimers(std::ostream& os);
};

#endif
