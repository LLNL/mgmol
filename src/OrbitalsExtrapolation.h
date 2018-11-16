// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#ifndef ORBITALSEXTRAPOLATION_H
#define ORBITALSEXTRAPOLATION_H

#include "Control.h"
class LocGridOrbitals;
class LocalizationRegions;
class ClusterOrbitals;
class MasksSet;
class HDFrestart;
class ProjectedMatricesInterface;

class OrbitalsExtrapolation
{
public:
    OrbitalsExtrapolation() { orbitals_minus1_ = 0; }

    virtual ~OrbitalsExtrapolation();

    virtual bool extrapolatedH() const { return false; }

    virtual void clearOldOrbitals();
    bool getRestartData(LocGridOrbitals& orbitals);
    virtual void setupPreviousOrbitals(LocGridOrbitals** orbitals,
        ProjectedMatricesInterface* proj_matrices, LocalizationRegions* lrs,
        ClusterOrbitals* local_cluster, MasksSet* currentMasks,
        MasksSet* corrtMasks, HDFrestart& h5f_file);

    virtual void extrapolate_orbitals(
        LocGridOrbitals** orbitals, LocGridOrbitals* new_orbitals)
        = 0;

    virtual short getNumAuxiliaryOrbitals() { return 0; }
    virtual short getNumOrbitalExtrapolations() { return 0; }

protected:
    LocGridOrbitals* orbitals_minus1_;
};

#endif
