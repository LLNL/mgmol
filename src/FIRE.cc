// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#include "FIRE.h"
#include "Control.h"
#include "Mesh.h"
#include "MasksSet.h"
#include "Electrostatic.h"
#include "KBPsiMatrix.h"
#include "MGmol.h"


FIRE::FIRE(LocGridOrbitals** orbitals,
           Ions& ions,
           Rho& rho,
           ConstraintSet& constraints,
           LocalizationRegions& lrs,
           MasksSet& masks,
           Electrostatic& electrostat,
           const double dt,
           MGmol& strategy):
              IonicAlgorithm(orbitals,
                             ions,
                             rho,
                             constraints,
                             lrs,
                             masks,
                             strategy),
    orbitals_(orbitals),
    ions_(ions),
    rho_(rho),
    lrs_(lrs),
    masks_(masks),
    electrostat_(electrostat),
    mgmol_strategy_(strategy)
{
    stepper_
        = new FIRE_IonicStepper(dt, atmove_, tau0_, taup_, fion_, 
                                 pmass_);
    
    registerStepper(stepper_);
}

FIRE::~FIRE()
{
}
