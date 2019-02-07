// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_RHO_H
#define MGMOL_RHO_H

#include <vector>

#include "Timer.h"
#include "DistMatrix.h"
#include "Control.h"
#include "global.h"

class HDFrestart;
class ProjectedMatricesInterface;

template <class T>
class Rho
{

    short nspin_;
    short myspin_;

    int np_;
#ifndef USE_DIS_MAT
    dist_matrix::DistMatrix<DISTMATDTYPE>& dm_;
#endif

    OrbitalsType orbitals_type_;

    std::vector<std::vector<int>> orbitals_indexes_;

    int iterative_index_;

    int verbosity_level_;

    // parameters for "blocking" in loops over functions and space
    int block_functions_;
    int block_space_;

    static Timer update_tm_;
    static Timer compute_tm_;
    static Timer compute_blas_tm_;

    double computeTotalCharge();
    void computeRhoSubdomain(const int iloc_init, const int iloc_end,
        const T& orbitals);
    void computeRhoSubdomain(const int iloc_init, const int iloc_end,
        const T& orbitals, const std::vector<DISTMATDTYPE>& occ);
    void computeRhoSubdomainOffDiagBlock(const int iloc_init,
        const int iloc_end,
        const std::vector<const T*>& vorbitals,
        const ProjectedMatricesInterface* const);
    void computeRhoSubdomainUsingBlas3(
        const int iloc_init, const int iloc_end, const T& orbitals);

    void accumulateCharge(const double alpha, const short ix_max,
        const ORBDTYPE* const psii, const ORBDTYPE* const psij,
        RHODTYPE* const plrho);
    int setupSubdomainData(const int iloc,
        const std::vector<const T*>& vorbitals,
        const ProjectedMatricesInterface* const projmatrices,
        std::vector<MATDTYPE>& melements,
        std::vector<std::vector<const ORBDTYPE*>>& mpsi);

    void computeRho(T& orbitals);
    void computeRho(
        T& orbitals, ProjectedMatricesInterface& proj_matrices);

public:
    // electronic density on grid
    std::vector<std::vector<RHODTYPE>> rho_;
    std::vector<std::vector<RHODTYPE>> rho_minus1_;

    Rho();
    ~Rho(){};

    void rescaleTotalCharge();
    void setup(const OrbitalsType orbitals_type,
               const std::vector<std::vector<int>>&);
    void setVerbosityLevel(const int vlevel) { verbosity_level_ = vlevel; }

    void update(T& current_orbitals);

    // compute rho using density matrix specified in arguments
    void computeRho(T& orbitals1, T& orbitals2,
        const dist_matrix::DistMatrix<DISTMATDTYPE>& dm11,
        const dist_matrix::DistMatrix<DISTMATDTYPE>& dm12,
        const dist_matrix::DistMatrix<DISTMATDTYPE>& dm21,
        const dist_matrix::DistMatrix<DISTMATDTYPE>& dm22);

    // compute rho using density matrix specified in arguments
    void computeRho(T& orbitals,
        const dist_matrix::DistMatrix<DISTMATDTYPE>& dm);

    void init(const RHODTYPE* const rhoc);
    void initUniform();
    int getIterativeIndex() const { return iterative_index_; }
    int readRestart(HDFrestart& file);
    void extrapolate();
    void axpyRhoc(const double alpha, RHODTYPE* rhoc);

    template <typename T2>
    double dotWithRho(const T2* const func) const;

    void gatherSpin();

    void setupBlockSizes(const int block_functions, const int block_space)
    {
        block_functions_ = block_functions;
        block_space_     = block_space;
    }

    static void printTimers(std::ostream& os);
};

#endif
