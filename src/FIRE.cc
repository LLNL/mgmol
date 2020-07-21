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

template <class T>
FIRE<T>::FIRE(T** orbitals, Ions& ions, Rho<T>& rho, ConstraintSet& constraints,
    std::shared_ptr<LocalizationRegions> lrs, MasksSet& masks,
    Electrostatic& electrostat, const double dt, MGmol<T>& strategy)
    : IonicAlgorithm<T>(orbitals, ions, rho, constraints, lrs, masks, strategy),
      orbitals_(orbitals),
      ions_(ions),
      rho_(rho),
      lrs_(lrs),
      masks_(masks),
      electrostat_(electrostat),
      mgmol_strategy_(strategy)
{
    stepper_ = new FIRE_IonicStepper(dt, IonicAlgorithm<T>::atmove_,
        IonicAlgorithm<T>::tau0_, IonicAlgorithm<T>::taup_,
        IonicAlgorithm<T>::fion_);

    IonicAlgorithm<T>::registerStepper(stepper_);
}

template class FIRE<LocGridOrbitals>;
template class FIRE<ExtendedGridOrbitals>;
