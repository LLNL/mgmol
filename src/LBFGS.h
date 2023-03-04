// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_LBFGS_H
#define MGMOL_LBFGS_H

#include "ClusterOrbitals.h"
#include "ConstraintSet.h"
#include "Energy.h"
#include "GridFunc.h"
#include "IonicAlgorithm.h"
#include "Ions.h"
#include "LBFGS_IonicStepper.h"
#include "LocalizationRegions.h"
#include "Rho.h"

class MasksSet;
class Electrostatic;
template <class OrbitalsType>
class MGmol;
class KBPsiMatrixInterface;

template <class OrbitalsType>
class LBFGS : public IonicAlgorithm<OrbitalsType>
{
private:
    OrbitalsType** orbitals_;
    Ions& ions_;
    Rho<OrbitalsType>& rho_;
    LBFGS_IonicStepper* stepper_;
    std::shared_ptr<LocalizationRegions> lrs_;
    ClusterOrbitals* local_cluster_;
    MasksSet& masks_;
    MasksSet& corrmasks_;
    Electrostatic& electrostat_;

    std::shared_ptr<LocalizationRegions> ref_lrs_;
    std::shared_ptr<MasksSet> ref_masks_;
    std::shared_ptr<MasksSet> ref_corrmasks_;
    std::shared_ptr<OrbitalsType> ref_orbitals_;
    pb::GridFunc<POTDTYPE>* vh_init_;

    double etot_i_[3];

    MGmol<OrbitalsType>& mgmol_strategy_;

    void updateRefMasks();
    void setup(const double dt);

public:
    LBFGS(OrbitalsType** orbitals, Ions& ions, Rho<OrbitalsType>& rho,
        ConstraintSet& constraints, std::shared_ptr<LocalizationRegions> lrs,
        ClusterOrbitals* local_cluster, MasksSet& masks, MasksSet& corrmasks,
        Electrostatic& electrostat, const double dt, MGmol<OrbitalsType>&);

    ~LBFGS() override;

    void reset(const double dt);

    int quenchElectrons(const int itmax, double& etot) override;
    void updateRefs();
    void setQuenchTol() const;
    void updatePotAndMasks();
    bool lbfgsLastStepNotAccepted() const;
};

#endif
