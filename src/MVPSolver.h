// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef MGMOL_MVPSOLVER_H
#define MGMOL_MVPSOLVER_H

#include "DistMatrix.h"
#include "Energy.h"
#include "MGmol.h"
#include "Rho.h"
#include "Timer.h"

class Ions;
class Electrostatic;
class ProjectedMatrices2N;
class ProjectedMatrices;

template <class T>
class MVPSolver
{
private:
    MPI_Comm comm_;
    std::ostream& os_;

    short n_inner_steps_;

    bool use_old_dm_;
    Ions& ions_;

    Rho<T>* rho_;
    Energy<T>* energy_;
    Electrostatic* electrostat_;

    int history_length_;
    std::vector<double> eks_history_;

    MGmol<T>* mgmol_strategy_;

    double de_old_;
    double de_;

    int numst_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* work_;
    ProjectedMatrices* proj_mat_work_;

    static Timer solve_tm_;
    static Timer target_tm_;

    double evaluateDerivative(dist_matrix::DistMatrix<DISTMATDTYPE>& dm2Ninit,
        dist_matrix::DistMatrix<DISTMATDTYPE>& delta_dm, const double ts0);
    void buildTarget_MVP(dist_matrix::DistMatrix<DISTMATDTYPE>& h11,
        dist_matrix::DistMatrix<DISTMATDTYPE>& s11,
        dist_matrix::DistMatrix<DISTMATDTYPE>& target);

public:
    MVPSolver(MPI_Comm comm, std::ostream& os, Ions& ions, Rho<T>* rho,
        Energy<T>* energy, Electrostatic* electrostat, MGmol<T>* mgmol_strategy,
        const int numst, const double kbT, const int nel,
        const std::vector<std::vector<int>>& global_indexes,
        const short n_inner_steps, const bool use_old_dm);
    ~MVPSolver();

    int solve(T& orbitals);
    void printTimers(std::ostream& os);
};

#endif
