// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_GRASSMANCG_H
#define MGMOL_GRASSMANCG_H

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrix.h"
#endif

#include "GrassmanLineMinimization.h"
#include "Hamiltonian.h"
#include "Potentials.h"
#include "ProjectedMatricesInterface.h"

#include <iostream>

template <class T>
class MGmol;

template <class T>
class GrassmanCG : public GrassmanLineMinimization<T>
{

public:
    GrassmanCG(Hamiltonian<T>* hamiltonian,
        ProjectedMatricesInterface* proj_matrices, MGmol<T>* mgmol_strategy,
        Ions& ions, std::ostream& os)
        : GrassmanLineMinimization<T>(
            hamiltonian, proj_matrices, mgmol_strategy, ions, os)
    {
    }

    void conjugate() override;
    double computeStepSize(T& orbitals) override;
#ifdef MGMOL_USE_SCALAPACK
    void computeOrbitalsProdWithH(
        T& orbitals1, T& orbitals2, dist_matrix::DistMatrix<DISTMATDTYPE>& mat);
    void computeOrbitalsProdWithH(
        T& orbitals, dist_matrix::DistMatrix<DISTMATDTYPE>& mat);
#endif
    void parallelTransportUpdate(const double lambda, T& orbitals) override;
};
#endif
