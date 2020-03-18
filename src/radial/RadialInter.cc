// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "RadialInter.h"

// Gregory-Newton cubic interpolation using forward differences
double RadialInter::cubint(const double r, const int j) const
{
    assert(y_.size() > 0);
    assert(j < (int)y_.size());
    assert(invdr_ > 0.);

    const std::vector<double>& yj = y_[j];
    double d0                     = r * invdr_;
    if (d0 < 1.)
    {
        return (1. - d0) * yj[0] + d0 * yj[1];
    }

    int ic = (int)d0;
    ic     = (ic > 0) ? ic : 1;
    if (ic + 2 >= (int)yj.size()) return 0.;

    d0 -= (double)(ic);
    double d1 = (d0 - 1.) * 0.5;
    double d2 = (d0 - 2.) / 3.;

    assert(yj[ic] < 10000.);
    assert(yj[ic] > -10000.);

    double f0 = yj[ic];
    double g0 = yj[ic] - yj[ic - 1];
    double g1 = yj[ic + 1] - yj[ic];
    double g2 = yj[ic + 2] - yj[ic + 1];
    double h1 = g1 - g0;
    double h2 = g2 - g1;
    double i2 = h2 - h1;

    return f0 + d0 * (g1 + d1 * (h2 + d2 * i2));
}

// linear interpolation
double RadialInter::linint(const double r, const int j) const
{
    assert(j < (int)y_.size());
    assert(invdr_ > 0.);

    const std::vector<double>& yj = y_[j];
    double b                      = r * invdr_;
    int ic                        = (int)b;
    if (ic >= (int)yj.size()) return 0.;

    b -= (double)(ic);
    double a = 1. - b;

    return a * yj[ic] + b * yj[ic + 1];
}
