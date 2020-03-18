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
template <class T>
class MGmol;
class KBPsiMatrixInterface;

template <class T>
class LBFGS : public IonicAlgorithm<T>
{
private:
    T** orbitals_;
    Ions& ions_;
    Rho<T>& rho_;
    LBFGS_IonicStepper* stepper_;
    LocalizationRegions& lrs_;
    ClusterOrbitals* local_cluster_;
    MasksSet& masks_;
    MasksSet& corrmasks_;
    Electrostatic& electrostat_;

    LocalizationRegions ref_lrs_;
    MasksSet* ref_masks_;
    MasksSet* ref_corrmasks_;
    T* ref_orbitals_;
    pb::GridFunc<POTDTYPE>* vh_init_;

    double etot_i_[3];

    MGmol<T>& mgmol_strategy_;

    void updateRefMasks();

public:
    LBFGS(T** orbitals, Ions& ions, Rho<T>& rho, ConstraintSet& constraints,
        LocalizationRegions& lrs, ClusterOrbitals* local_cluster,
        MasksSet& masks, MasksSet& corrmasks, Electrostatic& electrostat,
        const double dt, MGmol<T>&);

    ~LBFGS() override;

    int quenchElectrons(const int itmax, double& etot) override;
    void updateRefs();
    void setQuenchTol() const;
    void updatePotAndMasks();
    bool lbfgsLastStepNotAccepted() const;
};

#endif
