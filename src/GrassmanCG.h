// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef _GRASSMANCG_H_
#define _GRASSMANCG_H_

#include "Potentials.h"
#include "ProjectedMatricesInterface.h"
#include "GrassmanLineMinimization.h"
#include "SparseDistMatrix.h"
#include "DistMatrix.h"

#include <iostream>

class Energy;
class MGmol;

class GrassmanCG : public GrassmanLineMinimization
{      

public:
    GrassmanCG(Hamiltonian* hamiltonian,
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
    void computeOrbitalsProdWithH(LocGridOrbitals& orbitals1, LocGridOrbitals& orbitals2, dist_matrix::DistMatrix<DISTMATDTYPE>& mat);
    void computeOrbitalsProdWithH(LocGridOrbitals& orbitals, dist_matrix::DistMatrix<DISTMATDTYPE>& mat);  
    void parallelTransportUpdate(const double lambda, LocGridOrbitals& orbitals);
};
#endif  
