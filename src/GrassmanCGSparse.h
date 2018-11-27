// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_GRASSMANCGSPARSE_H
#define MGMOL_GRASSMANCGSPARSE_H

#include "GrassmanLineMinimization.h"
#include "Potentials.h"
#include "ProjectedMatricesInterface.h"
#include "VariableSizeMatrix.h"

#include <iostream>

class MGmol;

template <class T>
class GrassmanCGSparse : public GrassmanLineMinimization<T>
{

public:
    GrassmanCGSparse(Hamiltonian<T>* hamiltonian,
        ProjectedMatricesInterface* proj_matrices, MGmol* mgmol_strategy,
        Ions& ions, std::ostream& os)
        : GrassmanLineMinimization<T>(
              hamiltonian, proj_matrices, mgmol_strategy, ions, os)
    {
    }

    void conjugate();
    double computeStepSize(T& orbitals);
    void computeOrbitalsProdWithH(T& orbitals1,
        T& orbitals2, VariableSizeMatrix<sparserow>& mat,
        const bool consolidate);
    void computeOrbitalsProdWithH(T& orbitals,
        VariableSizeMatrix<sparserow>& mat, const bool consolidate);
    void parallelTransportUpdate(
        const double lambda, T& orbitals);
};
#endif
