// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
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

#ifndef FIRE_IonicStepper_H
#define FIRE_IonicStepper_H

#include "IonicStepper.h"
#include <cassert>

class FIRE_IonicStepper : public IonicStepper
{
private:
    std::vector<double>& vel_;

    std::vector<double>& fion_;

    double dt_;

    ///
    /// Inverse of largest mass in system
    ///
    double invmass_;

    ///
    /// FIRE algorithm parameters
    ///
    int nmin_;
    int npp_;

    double alpha_init_;
    double alpha_;
    double falpha_;
    double dtmax_;
    double finc_;
    double fdec_;

public:
    FIRE_IonicStepper(const double dt, const std::vector<short>& atmove,
        std::vector<double>& tau0, std::vector<double>& taup,
        std::vector<double>& vel, std::vector<double>& fion,
        std::vector<double>& masses);

    int run() override;
    double etol(void) const override;
    int write_hdf5(HDFrestart&) override;
    int init(HDFrestart&) override;
};

#endif
