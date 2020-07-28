// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "DistMatrix.h"
#include "ProjectedMatrices2N.h"
#include "Timer.h"

class Ions;
template <class T>
class Hamiltonian;
template <class T>
class Rho;
template <class T>
class Energy;
template <class T>
class MGmol;
class Electrostatic;

template <class T>
class DavidsonSolver
{
private:
    MPI_Comm comm_;
    std::ostream& os_;

    Ions& ions_;

    Hamiltonian<T>* hamiltonian_;
    Rho<T>* rho_;
    Energy<T>* energy_;
    Electrostatic* electrostat_;

    int history_length_;
    std::vector<double> eks_history_;

    MGmol<T>* mgmol_strategy_;

    double de_old_;
    double de_;

    int numst_;
    std::unique_ptr<dist_matrix::DistMatrix<DISTMATDTYPE>> work2N_;
    std::unique_ptr<ProjectedMatrices2N<dist_matrix::DistMatrix<DISTMATDTYPE>>>
        proj_mat2N_;

    static Timer solve_tm_;
    static Timer target_tm_;

    // void swapColumnsVect(dist_matrix::DistMatrix<DISTMATDTYPE>& evect,
    //    const dist_matrix::DistMatrix<DISTMATDTYPE>& hb2N,
    //    const std::vector<DISTMATDTYPE>& eval,
    //    dist_matrix::DistMatrix<DISTMATDTYPE>& work2N);
    int checkConvergence(const double e0, const int it, const double tol);
    double evaluateDerivative(dist_matrix::DistMatrix<DISTMATDTYPE>& dm2Ninit,
        dist_matrix::DistMatrix<DISTMATDTYPE>& delta_dm, const double ts0);
    void buildTarget2N_MVP(dist_matrix::DistMatrix<DISTMATDTYPE>& h11,
        dist_matrix::DistMatrix<DISTMATDTYPE>& h12,
        dist_matrix::DistMatrix<DISTMATDTYPE>& h21,
        dist_matrix::DistMatrix<DISTMATDTYPE>& h22,
        dist_matrix::DistMatrix<DISTMATDTYPE>& s11,
        dist_matrix::DistMatrix<DISTMATDTYPE>& s22,
        dist_matrix::DistMatrix<DISTMATDTYPE>& target);
    // void buildTarget2N_new(dist_matrix::DistMatrix<DISTMATDTYPE>& h11,
    //    dist_matrix::DistMatrix<DISTMATDTYPE>& h12,
    //    dist_matrix::DistMatrix<DISTMATDTYPE>& h21,
    //    dist_matrix::DistMatrix<DISTMATDTYPE>& h22,
    //    dist_matrix::DistMatrix<DISTMATDTYPE>& s11,
    //    dist_matrix::DistMatrix<DISTMATDTYPE>& s22,
    //    const std::vector<DISTMATDTYPE>& occ,
    //    const std::vector<DISTMATDTYPE>& auxenergies, const double kbT,
    //    const double eta, dist_matrix::DistMatrix<DISTMATDTYPE>& target);

public:
    DavidsonSolver(MPI_Comm comm, std::ostream& os, Ions& ions,
        Hamiltonian<T>* hamiltonian, Rho<T>* rho, Energy<T>* energy,
        Electrostatic* electrostat, MGmol<T>* mgmol_strategy, const int numst,
        const double kbT, const int nel,
        const std::vector<std::vector<int>>& global_indexes);
    ~DavidsonSolver();

    int solve(T& orbitals, T& work_orbitals);

    static void printTimers(std::ostream&);
};
