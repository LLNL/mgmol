// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "DielectricControl.h"
#include "ProjectedMatrices2N.h"
#include "Timer.h"

class Ions;
template <class OrbitalsType>
class Hamiltonian;
template <class OrbitalsType>
class Rho;
template <class OrbitalsType>
class Energy;
template <class OrbitalsType>
class MGmol;
class Electrostatic;

template <class OrbitalsType, class MatrixType>
class DavidsonSolver
{
private:
    std::ostream& os_;

    Ions& ions_;

    Hamiltonian<OrbitalsType>* hamiltonian_;
    Rho<OrbitalsType>* rho_;
    Energy<OrbitalsType>* energy_;
    Electrostatic* electrostat_;

    int history_length_;
    std::vector<double> eks_history_;

    MGmol<OrbitalsType>* mgmol_strategy_;

    double de_old_;
    double de_;
    double mixing_;

    int numst_;
    std::unique_ptr<MatrixType> work2N_;
    std::unique_ptr<ProjectedMatrices2N<MatrixType>> proj_mat2N_;

    DielectricControl diel_control_;

    static Timer solve_tm_;
    static Timer target_tm_;

    // void swapColumnsVect(MatrixType& evect,
    //    const MatrixType& hb2N,
    //    const std::vector<DISTMATDTYPE>& eval,
    //    MatrixType& work2N);
    int checkConvergence(const double e0, const int it, const double tol);
    double evaluateDerivative(
        MatrixType& dm2Ninit, MatrixType& delta_dm, const double ts0);
    void buildTarget2N_MVP(MatrixType& h11, MatrixType& h12, MatrixType& h21,
        MatrixType& h22, MatrixType& s11, MatrixType& s22, MatrixType& target);
    // void buildTarget2N_new(MatrixType& h11,
    //    MatrixType& h12,
    //    MatrixType& h21,
    //    MatrixType& h22,
    //    MatrixType& s11,
    //    MatrixType& s22,
    //    const std::vector<DISTMATDTYPE>& occ,
    //    const std::vector<DISTMATDTYPE>& auxenergies, const double kbT,
    //    const double eta, MatrixType& target);

public:
    DavidsonSolver(std::ostream& os, Ions& ions,
        Hamiltonian<OrbitalsType>* hamiltonian, Rho<OrbitalsType>* rho,
        Energy<OrbitalsType>* energy, Electrostatic* electrostat,
        MGmol<OrbitalsType>* mgmol_strategy,
        const std::vector<std::vector<int>>& global_indexes,
        const double mixing, const bool with_spin);
    ~DavidsonSolver();

    int solve(OrbitalsType& orbitals, OrbitalsType& work_orbitals);

    static void printTimers(std::ostream&);
};
