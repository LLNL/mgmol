// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef _GRASSMANCGSPARSE_H_
#define _GRASSMANCGSPARSE_H_

#include "Potentials.h"
#include "ProjectedMatricesInterface.h"
#include "GrassmanLineMinimization.h"
#include "VariableSizeMatrix.h"

#include <iostream>

class MGmol;

class GrassmanCGSparse : public GrassmanLineMinimization
{      
    
public:
    GrassmanCGSparse(Hamiltonian* hamiltonian,
         ProjectedMatricesInterface* proj_matrices,
         MGmol* mgmol_strategy,
         Ions& ions,
         std::ostream& os)
        : GrassmanLineMinimization(hamiltonian,
          proj_matrices,
          mgmol_strategy,
          ions,
          os)
    {
    }

    void conjugate();
    double computeStepSize(LocGridOrbitals& orbitals);                         
    void computeOrbitalsProdWithH(LocGridOrbitals& orbitals1, LocGridOrbitals& orbitals2, VariableSizeMatrix<sparserow>& mat, const bool consolidate);
    void computeOrbitalsProdWithH(LocGridOrbitals& orbitals, VariableSizeMatrix<sparserow>& mat, const bool consolidate);  
    void parallelTransportUpdate(const double lambda, LocGridOrbitals& orbitals);
};
#endif  
