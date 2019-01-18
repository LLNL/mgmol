// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ORBITALSEXTRAPOLATIONORDER3_H
#define MGMOL_ORBITALSEXTRAPOLATIONORDER3_H

#include "OrbitalsExtrapolation.h"

template <class T>
class OrbitalsExtrapolationOrder3 : public OrbitalsExtrapolation<T>
{
private:
    T* initial_orbitals_minus2_;
    T* orbitals_minus1_;
    T* orbitals_minus2_;

public:
    OrbitalsExtrapolationOrder3()
        : initial_orbitals_minus2_(0),
          orbitals_minus1_(0),
          orbitals_minus2_(0)
    {
    };

    ~OrbitalsExtrapolationOrder3()
    {
        if (orbitals_minus2_ != 0)
        {
            delete orbitals_minus2_;
            orbitals_minus2_ = 0;
        }
        if (initial_orbitals_minus2_ != 0)
        {
            delete initial_orbitals_minus2_;
            initial_orbitals_minus2_ = 0;
        }
    }

    void extrapolate_orbitals(
        T** orbitals, T* new_orbitals);

    void clearOldOrbitals()
    {
        OrbitalsExtrapolation<T>::clearOldOrbitals();

        if (orbitals_minus2_ != 0)
        {
            delete orbitals_minus2_;
            orbitals_minus2_ = 0;
        }
    }
};

#endif
