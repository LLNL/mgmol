// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#ifndef FIRE_H
#define FIRE_H

#include "Energy.h"
#include "FIRE_IonicStepper.h"
#include "IonicAlgorithm.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "Rho.h"

class MasksSet;
class Energy;
class Electrostatic;
class MGmol;
class KBPsiMatrixInterface;
class ConstraintSet;

class FIRE : public IonicAlgorithm
{
private:
    LocGridOrbitals** orbitals_;
    const Ions& ions_;
    const Rho& rho_;
    FIRE_IonicStepper* stepper_;
    const LocalizationRegions& lrs_;
    const MasksSet& masks_;
    const Electrostatic& electrostat_;

    const MGmol& mgmol_strategy_;

public:
    FIRE(LocGridOrbitals** orbitals, Ions& ions, Rho& rho,
        ConstraintSet& constraints, LocalizationRegions& lrs, MasksSet& masks,
        Electrostatic& electrostat, const double dt, MGmol&);

    ~FIRE();
};

#endif
