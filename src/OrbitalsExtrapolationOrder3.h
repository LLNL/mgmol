// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ORBITALSEXTRAPOLATIONORDER3_H
#define MGMOL_ORBITALSEXTRAPOLATIONORDER3_H

#include "OrbitalsExtrapolation.h"

template <class OrbitalsType>
class OrbitalsExtrapolationOrder3 : public OrbitalsExtrapolation<OrbitalsType>
{
private:
    OrbitalsType* initial_orbitals_minus2_;
    OrbitalsType* orbitals_minus1_;
    OrbitalsType* orbitals_minus2_;

public:
    OrbitalsExtrapolationOrder3()
        : initial_orbitals_minus2_(nullptr),
          orbitals_minus1_(nullptr),
          orbitals_minus2_(nullptr) {};

    ~OrbitalsExtrapolationOrder3() override
    {
        if (orbitals_minus2_ != nullptr)
        {
            delete orbitals_minus2_;
            orbitals_minus2_ = nullptr;
        }
        if (initial_orbitals_minus2_ != nullptr)
        {
            delete initial_orbitals_minus2_;
            initial_orbitals_minus2_ = nullptr;
        }
    }

    void extrapolate_orbitals(
        OrbitalsType** orbitals, OrbitalsType* new_orbitals) override;

    void clearOldOrbitals() override
    {
        OrbitalsExtrapolation<OrbitalsType>::clearOldOrbitals();

        if (orbitals_minus2_ != nullptr)
        {
            delete orbitals_minus2_;
            orbitals_minus2_ = nullptr;
        }
    }
};

#endif
