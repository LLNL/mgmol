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

#include "Energy.h"
#include "FIRE_IonicStepper.h"
#include "IonicAlgorithm.h"
#include "Ions.h"
#include "LocalizationRegions.h"
#include "MGmol.h"
#include "Rho.h"

class MasksSet;
class Electrostatic;
class KBPsiMatrixInterface;
class ConstraintSet;

template <class T>
class FIRE : public IonicAlgorithm<T>
{
private:
    T** orbitals_;
    const Ions& ions_;
    const Rho<T>& rho_;
    FIRE_IonicStepper* stepper_;
    const std::shared_ptr<LocalizationRegions> lrs_;
    const MasksSet& masks_;
    const Electrostatic& electrostat_;

    const MGmol<T>& mgmol_strategy_;

public:
    FIRE(T** orbitals, Ions& ions, Rho<T>& rho, ConstraintSet& constraints,
        std::shared_ptr<LocalizationRegions> lrs, MasksSet& masks,
        Electrostatic& electrostat, const double dt, MGmol<T>&);

    ~FIRE() override{};
};

#endif
