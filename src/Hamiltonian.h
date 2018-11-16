// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_HAMILTONIAN_H
#define MGMOL_HAMILTONIAN_H

#include "LapFactory.h"
#include "LocGridOrbitals.h"
#include "SparseDistMatrix.h"
#include "Timer.h"
#include "VariableSizeMatrix.h"

class LocGridOrbitals;
class Potentials;

class Hamiltonian
{

    pb::Lap<ORBDTYPE>* lapOper_;
    Potentials* pot_;
    LocGridOrbitals* hlphi_;
    int itindex_;

    static Timer apply_Hloc_tm_;

    void applyLocal(const int istate, const int nstates, LocGridOrbitals& phi,
        LocGridOrbitals& hphi);

public:
    static Timer apply_Hloc_tm() { return apply_Hloc_tm_; }

    Hamiltonian();
    ~Hamiltonian();

    void setup(const pb::Grid& myGrid, const int lap_type);

    Potentials& potential() { return *pot_; }
    void setHlOutdated() { itindex_ = -1; }
    pb::Lap<ORBDTYPE>* lapOper() { return lapOper_; }

    const LocGridOrbitals& applyLocal(
        LocGridOrbitals& phi, const bool force = false);

    void addHlocal2matrix(LocGridOrbitals& orbitals1,
        LocGridOrbitals& orbitals2,
        dist_matrix::SparseDistMatrix<DISTMATDTYPE>& mat,
        const bool force = false);
    void addHlocal2matrix(LocGridOrbitals& orbitals1,
        LocGridOrbitals& orbitals2, VariableSizeMatrix<sparserow>& mat,
        const bool force = false);
    void addHlocalij(LocGridOrbitals& orbitals1, LocGridOrbitals& orbitals2,
        ProjectedMatricesInterface*);
};

#endif
