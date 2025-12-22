// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_OrbitalsPreconditioning_H
#define MGMOL_OrbitalsPreconditioning_H

#include "LocalizationRegions.h"
#include "MasksSet.h"

template <class OrbitalsType>
class OrbitalsPreconditioning
{
public:
    OrbitalsPreconditioning() {};

    virtual ~OrbitalsPreconditioning() {};

    virtual void setup(OrbitalsType& orbitals, MasksSet*,
        const std::shared_ptr<LocalizationRegions>&)
        = 0;

    virtual void precond(OrbitalsType& orbitals) = 0;
};

#endif
