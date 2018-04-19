// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef LBFGS_H
#define LBFGS_H

#include "IonicAlgorithm.h"
#include "LBFGS_IonicStepper.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "Rho.h"
#include "ConstraintSet.h"
#include "Energy.h"
#include "LocalizationRegions.h"
#include "GridFunc.h"
#include "ClusterOrbitals.h"

class MasksSet;
class Electrostatic;
class MGmol;
class KBPsiMatrixInterface;

class LBFGS : public IonicAlgorithm
{
private:
    LocGridOrbitals** orbitals_;
    Ions& ions_;
    Rho& rho_;
    LBFGS_IonicStepper* stepper_;
    LocalizationRegions& lrs_;
    ClusterOrbitals* local_cluster_;
    MasksSet& masks_;
    MasksSet& corrmasks_;
    Electrostatic& electrostat_;
    
    LocalizationRegions ref_lrs_;
    MasksSet* ref_masks_;
    MasksSet* ref_corrmasks_;
    LocGridOrbitals* ref_orbitals_;
    pb::GridFunc<POTDTYPE>* vh_init_;

    double etot_i_[3];
    
    MGmol& mgmol_strategy_;

public:
    LBFGS(LocGridOrbitals** orbitals,
          Ions& ions,
          Rho& rho,
          ConstraintSet& constraints,
          LocalizationRegions& lrs,
          ClusterOrbitals* local_cluster,
          MasksSet& masks,
          MasksSet& corrmasks,
          Electrostatic& electrostat,
          const double dt,
          MGmol&);
    
    ~LBFGS();
    
    int quenchElectrons(const int itmax, double& etot);
    void updateRefs();
    void setQuenchTol()const;
    void updatePotAndMasks();
    bool lbfgsLastStepNotAccepted()const;
};

#endif
