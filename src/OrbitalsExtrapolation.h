// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ORBITALSEXTRAPOLATION_H
#define MGMOL_ORBITALSEXTRAPOLATION_H

#include <memory>

class LocalizationRegions;
class ClusterOrbitals;
class MasksSet;
class HDFrestart;
class ProjectedMatricesInterface;

template <class OrbitalsType>
class OrbitalsExtrapolation
{
public:
    OrbitalsExtrapolation() : orbitals_minus1_(nullptr) { }

    virtual ~OrbitalsExtrapolation();

    virtual void clearOldOrbitals();
    bool getRestartData(OrbitalsType& orbitals);
    virtual void setupPreviousOrbitals(OrbitalsType** orbitals,
        ProjectedMatricesInterface* proj_matrices,
        std::shared_ptr<LocalizationRegions> lrs,
        ClusterOrbitals* local_cluster, MasksSet* currentMasks,
        MasksSet* corrtMasks, HDFrestart& h5f_file);

    virtual void extrapolate_orbitals(
        OrbitalsType** orbitals, OrbitalsType* new_orbitals)
        = 0;

    virtual short getNumAuxiliaryOrbitals() { return 0; }
    virtual short getNumOrbitalExtrapolations() { return 0; }

protected:
    OrbitalsType* orbitals_minus1_ = nullptr;
};

#endif
