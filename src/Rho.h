// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef RHO_H
#define RHO_H

#include <vector>

// pb
#include "Timer.h"

#include "DistMatrix.h"
#include "LocGridOrbitals.h"

#include "global.h"

//class LocGridOrbitals;
class HDFrestart;
class ProjectedMatricesInterface;


class Rho{

    short nspin_;
    short myspin_;
    
    int np_;
#ifndef USE_DIS_MAT
    dist_matrix::DistMatrix<DISTMATDTYPE>& dm_;
#endif
    
    int orbitals_type_;
    
    std::vector<std::vector<int> > orbitals_indexes_;

    int iterative_index_;
    
    int verbosity_level_;

    //parameters for "blocking" in loops over functions and space
    int block_functions_;
    int block_space_;
    
    static Timer update_tm_;
    static Timer compute_tm_;

    double computeTotalCharge();
    void computeRhoSubdomain(const int iloc_init, 
                            const int iloc_end,
                            const LocGridOrbitals& orbitals);
    void computeRhoSubdomain(const int iloc_init, 
                             const int iloc_end,
                             const LocGridOrbitals& orbitals,
                             const std::vector<DISTMATDTYPE>& occ);
    void computeRhoSubdomainOffDiagBlock(const int iloc_init, 
                                        const int iloc_end,
                                        const std::vector<const LocGridOrbitals*>& vorbitals,
                                        const ProjectedMatricesInterface* const);

    void accumulateCharge(
        const double alpha,
        const short ix_max,
        const ORBDTYPE* const psii,
        const ORBDTYPE* const psij,
        RHODTYPE* const plrho);
    int setupSubdomainData(const int iloc,
                           const std::vector<const LocGridOrbitals*>& vorbitals,
                           const ProjectedMatricesInterface* const projmatrices,
                           std::vector<MATDTYPE>& melements,
                           std::vector<std::vector<const ORBDTYPE*> >& mpsi);

    void computeRho(LocGridOrbitals& orbitals);
    void computeRho(LocGridOrbitals& orbitals, ProjectedMatricesInterface& proj_matrices);

public:

    // electronic density on grid
    std::vector<std::vector<RHODTYPE> > rho_;
    std::vector<std::vector<RHODTYPE> > rho_minus1_;
    
    Rho();
    ~Rho();
   
    void rescaleTotalCharge();
    void setup(const int orbitals_type,
               const std::vector<std::vector<int> >&);
    void setVerbosityLevel(const int vlevel)
    {
        verbosity_level_=vlevel;
    }
    
    void update(LocGridOrbitals& current_orbitals);

    //compute rho using density matrix specified in arguments
    void computeRho(LocGridOrbitals& orbitals1,
                    LocGridOrbitals& orbitals2,
                    const dist_matrix::DistMatrix<DISTMATDTYPE>& dm11,
                    const dist_matrix::DistMatrix<DISTMATDTYPE>& dm12,
                    const dist_matrix::DistMatrix<DISTMATDTYPE>& dm21,
                    const dist_matrix::DistMatrix<DISTMATDTYPE>& dm22);

    //compute rho using density matrix specified in arguments
    void computeRho(LocGridOrbitals& orbitals,
                    const dist_matrix::DistMatrix<DISTMATDTYPE>& dm);

    void init(const RHODTYPE* const rhoc);
    void initUniform();
    int getIterativeIndex()const{ return iterative_index_; }
    int readRestart(HDFrestart& file);
    void extrapolate();
    void axpyRhoc(const double alpha, RHODTYPE *rhoc);
    
    template <typename T>
    double dotWithRho(const T* const func)const;
    
    void gatherSpin();
    
    void setupBlockSizes(const int block_functions, const int block_space)
    {
        block_functions_=block_functions;
        block_space_=block_space;
    }

    static void printTimers(std::ostream& os);
};

#endif
