// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_OrbitalsExtrapolationFACTORY_H
#define MGMOL_OrbitalsExtrapolationFACTORY_H

#include "OrbitalsExtrapolationOrder2.h"
#include "OrbitalsExtrapolationOrder3.h"
#include "SpreadPenalty.h"

template <class T>
class OrbitalsExtrapolationFactory
{
public:
    static OrbitalsExtrapolation<T>* create(
        const WFExtrapolationType type, SpreadPenalty<T>* spread_penalty)
    {
        OrbitalsExtrapolation<T>* orbitals_extrapol;
        switch (type)
        {
            case WFExtrapolationType::Order2:
                orbitals_extrapol = new OrbitalsExtrapolationOrder2<T>();
                break;
            case WFExtrapolationType::Order3:
                orbitals_extrapol = new OrbitalsExtrapolationOrder3<T>();
                break;
            default:
                (*MPIdata::serr)
                    << "OrbitalsExtrapolation* create() --- option invalid\n";
                exit(2);
        }
        return orbitals_extrapol;
    }
};

#endif
