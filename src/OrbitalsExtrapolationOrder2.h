// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ORBITALSEXTRAPOLATIONORDER2_H
#define MGMOL_ORBITALSEXTRAPOLATIONORDER2_H

#include "OrbitalsExtrapolation.h"

template <class OrbitalsType>
class OrbitalsExtrapolationOrder2 : public OrbitalsExtrapolation<OrbitalsType>
{
public:
    OrbitalsExtrapolationOrder2() { }

    void extrapolate_orbitals(
        OrbitalsType** orbitals, OrbitalsType* new_orbitals) override;
};

#endif
