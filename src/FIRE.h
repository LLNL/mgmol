// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_FIRE_H
#define MGMOL_FIRE_H

#include "ConstraintSet.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "FIRE_IonicStepper.h"
#include "IonicAlgorithm.h"
#include "Ions.h"
#include "LocalizationRegions.h"
#include "MGmol.h"
#include "MasksSet.h"
#include "Rho.h"

template <class OrbitalsType>
class FIRE : public IonicAlgorithm<OrbitalsType>
{
private:
    OrbitalsType** orbitals_;
    const Ions& ions_;
    const Rho<OrbitalsType>& rho_;
    FIRE_IonicStepper* stepper_;
    const std::shared_ptr<LocalizationRegions> lrs_;
    const MasksSet& masks_;
    const Electrostatic& electrostat_;

    const MGmol<OrbitalsType>& mgmol_strategy_;

public:
    FIRE(OrbitalsType** orbitals, Ions& ions, Rho<OrbitalsType>& rho,
        ConstraintSet& constraints, std::shared_ptr<LocalizationRegions> lrs,
        MasksSet& masks, Electrostatic& electrostat, const double dt,
        MGmol<OrbitalsType>&);

    ~FIRE() override {};
};

#endif
