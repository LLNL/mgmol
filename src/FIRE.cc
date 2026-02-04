// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "FIRE.h"
#include "Control.h"
#include "Electrostatic.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "MasksSet.h"
#include "Mesh.h"

template <class OrbitalsType>
FIRE<OrbitalsType>::FIRE(OrbitalsType** orbitals, Ions& ions,
    Rho<OrbitalsType>& rho, ConstraintSet& constraints,
    std::shared_ptr<LocalizationRegions> lrs, MasksSet& masks,
    Electrostatic& electrostat, const double dt, MGmol<OrbitalsType>& strategy)
    : IonicAlgorithm<OrbitalsType>(
        orbitals, ions, rho, constraints, lrs, masks, strategy),
      orbitals_(orbitals),
      ions_(ions),
      rho_(rho),
      lrs_(lrs),
      masks_(masks),
      electrostat_(electrostat),
      mgmol_strategy_(strategy)
{
    stepper_ = new FIRE_IonicStepper(dt, IonicAlgorithm<OrbitalsType>::atmove_,
        IonicAlgorithm<OrbitalsType>::tau0_,
        IonicAlgorithm<OrbitalsType>::taup_, ions.getVelocities(),
        IonicAlgorithm<OrbitalsType>::fion_,
        IonicAlgorithm<OrbitalsType>::pmass_);

    IonicAlgorithm<OrbitalsType>::registerStepper(stepper_);
}

template class FIRE<LocGridOrbitals<ORBDTYPE>>;
template class FIRE<ExtendedGridOrbitals<ORBDTYPE>>;
