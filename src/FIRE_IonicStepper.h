// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

////////////////////////////////////////////////////////////////////////////////
//
// FIRE algorithm
// according to Bitzek et al., Phys. rev. Lett. 97, 170201 (2006)
//
////////////////////////////////////////////////////////////////////////////////
// $Id:$

#ifndef FIRE_IonicStepper_H
#define FIRE_IonicStepper_H

#include "IonicStepper.h"
#include <cassert>

class FIRE_IonicStepper : public IonicStepper
{
private:
    // taum_ stores velocities if dt_<=0
    std::vector<double> taum_; // taum_[3*na_]

    const std::vector<double>& pmass_; // pmass[ia]
    std::vector<double>& fion_;

    int nmin_;
    int npp_;

    double alpha_init_;
    double alpha_;
    double falpha_;
    double dt_;
    double dtmax_;
    double finc_;
    double fdec_;

public:
    FIRE_IonicStepper(const double dt, const std::vector<short>& atmove,
        std::vector<double>& tau0, std::vector<double>& taup,
        std::vector<double>& fion, const std::vector<double>& pmass);

    int run();
    double etol(void) const;
    int write_hdf5(HDFrestart&);
    int init(HDFrestart&);
};

#endif
