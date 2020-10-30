// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_RHO_H
#define MGMOL_RHO_H

#include "Control.h"
#include "Timer.h"
#include "global.h"

#include <vector>

class HDFrestart;
class ProjectedMatricesInterface;

template <class OrbitalsType>
class Rho
{
    short nspin_;
    short myspin_;

    int np_;

    OrthoType orbitals_type_;

    std::vector<std::vector<int>> orbitals_indexes_;

    int iterative_index_;

    int verbosity_level_;

    // parameters for "blocking" in loops over functions and space
    // int block_functions_;
    // int block_space_;

    static Timer update_tm_;
    static Timer compute_tm_;
    static Timer compute_blas_tm_;
    // static Timer compute_offdiag_tm_;

    double computeTotalCharge();
    // void computeRhoSubdomain(const int iloc_init, const int iloc_end,
    //    const T& orbitals);
    void computeRhoSubdomain(const int iloc_init, const int iloc_end,
        const OrbitalsType& orbitals, const std::vector<double>& occ);
    // void computeRhoSubdomainOffDiagBlock(const int iloc_init,
    //    const int iloc_end,
    //    const std::vector<const T*>& vorbitals,
    //    const ProjectedMatricesInterface* const);
    void computeRhoSubdomainUsingBlas3(const int iloc_init, const int iloc_end,
        const OrbitalsType& orbitals1, const OrbitalsType& orbitals2);
    void computeRhoSubdomainUsingBlas3(
        const int iloc_init, const int iloc_end, const OrbitalsType& orbitals)
    {
        computeRhoSubdomainUsingBlas3(iloc_init, iloc_end, orbitals, orbitals);
    }

    void accumulateCharge(const double alpha, const short ix_max,
        const ORBDTYPE* const psii, const ORBDTYPE* const psij,
        RHODTYPE* const plrho);
    // int setupSubdomainData(const int iloc,
    //    const std::vector<const T*>& vorbitals,
    //    const ProjectedMatricesInterface* const projmatrices,
    //    std::vector<MATDTYPE>& melements,
    //    std::vector<std::vector<const ORBDTYPE*>>& mpsi);

    void computeRho(OrbitalsType& orbitals);
    void computeRho(
        OrbitalsType& orbitals, ProjectedMatricesInterface& proj_matrices);

    void gatherSpin();

public:
    // electronic density on grid
    std::vector<std::vector<RHODTYPE>> rho_;
    std::vector<std::vector<RHODTYPE>> rho_minus1_;

    Rho();
    ~Rho(){};

    void rescaleTotalCharge();
    void setup(
        const OrthoType orbitals_type, const std::vector<std::vector<int>>&);
    void setVerbosityLevel(const int vlevel) { verbosity_level_ = vlevel; }

    void update(OrbitalsType& current_orbitals);

    // compute rho using density matrix specified in arguments
    template <class MatrixType>
    void computeRho(OrbitalsType& orbitals1, OrbitalsType& orbitals2,
        const MatrixType& dm11, const MatrixType& dm12, const MatrixType& dm21,
        const MatrixType& dm22);

    // compute rho using density matrix specified in arguments
    template <class MatrixType>
    void computeRho(OrbitalsType& orbitals, const MatrixType& dm);

    void init(const RHODTYPE* const rhoc);
    void initUniform();
    int getIterativeIndex() const { return iterative_index_; }
    int readRestart(HDFrestart& file);
    void extrapolate();
    void axpyRhoc(const double alpha, RHODTYPE* rhoc);

    template <typename T2>
    double dotWithRho(const T2* const func) const;

    // void setupBlockSizes(const int block_functions, const int block_space)
    //{
    //    block_functions_ = block_functions;
    //    block_space_     = block_space;
    //}

    static void printTimers(std::ostream& os);
};

#endif
