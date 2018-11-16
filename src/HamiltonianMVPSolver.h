// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef HAMILTONIANMVP_SOLVER_H_
#define HAMILTONIANMVP_SOLVER_H_
// $Id:$
#include "DistMatrix.h"
// pb
#include "Timer.h"

class LocGridOrbitals;
class Ions;
class Rho;
class Energy;
class MGmol;
class Electrostatic;
class ProjectedMatrices2N;
class ProjectedMatrices;

template <class T1, class T2, class T3>
class HamiltonianMVPSolver
{

private:
    MPI_Comm comm_;
    std::ostream& os_;

    short n_inner_steps_;

    Ions& ions_;

    Rho* rho_;
    Energy* energy_;
    Electrostatic* electrostat_;
    MGmol* mgmol_strategy_;

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
    HamiltonianMVPSolver(MPI_Comm comm, std::ostream& os, Ions& ions, Rho* rho,
        Energy* energy, Electrostatic* electrostat, MGmol* mgmol_strategy,
        const int numst, const double kbT, const int nel,
        const std::vector<std::vector<int>>& global_indexes,
        const short n_inner_steps, const T1& hinit,
        const bool try_shorter_intervals = false);
    ~HamiltonianMVPSolver();
    int solve(LocGridOrbitals& orbitals);
    void reset();
    void printTimers(std::ostream& os);
};

#endif
