// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ORBITALSEXTRAPOLATION_H
#define MGMOL_ORBITALSEXTRAPOLATION_H

#define EXTRAPOLATE_H 1

#include "Hextrapolation.h"
#include "Control.h"
class LocalizationRegions;
class ClusterOrbitals;
class MasksSet;
class HDFrestart;
class ProjectedMatricesInterface;

template <class T>
class OrbitalsExtrapolation
{
public:
    OrbitalsExtrapolation()
        : orbitals_minus1_(0), hextrapol_(0)
    {
#if EXTRAPOLATE_H
        Control& ct = *(Control::instance());

        if (ct.OuterSolver() != OuterSolverType::ABPG &&
            ct.OuterSolver() != OuterSolverType::NLCG )
            hextrapol_ = new Hextrapolation(ct.numst);
#endif
    }

    virtual ~OrbitalsExtrapolation();

#if EXTRAPOLATE_H
    void initExtrapolationH(const dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
    {
        assert( hextrapol_ );

        hextrapol_->initExtrapolationH(mat);
    }
#endif

    virtual bool extrapolatedH() const { return false; }

    virtual void clearOldOrbitals();
    bool getRestartData(T& orbitals);
    virtual void setupPreviousOrbitals(T** orbitals,
        ProjectedMatricesInterface* proj_matrices, LocalizationRegions* lrs,
        ClusterOrbitals* local_cluster, MasksSet* currentMasks,
        MasksSet* corrtMasks, HDFrestart& h5f_file);

    virtual void extrapolate_orbitals(
        T** orbitals, T* new_orbitals)
        = 0;

    virtual short getNumAuxiliaryOrbitals() { return 0; }
    virtual short getNumOrbitalExtrapolations() { return 0; }

protected:
    T* orbitals_minus1_;

#if EXTRAPOLATE_H
    Hextrapolation* hextrapol_;
#endif
};

#endif
