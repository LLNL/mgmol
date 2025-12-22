// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "RadialSpline.h"

// initialize y2_
void RadialSpline::spline(const double yp1, const double ypn)
{
    const int n = (int)x_.size();
    vector<double> u(n - 1);

    y2_.resize(x_.size());

    if (yp1 >= 1.e30)
    {
        y2_[0] = 0.0;
        u[0]   = 0.0;
    }
    else
    {
        y2_[0] = -0.5;
        assert(x_[1] - x_[0] > 0.0);
        u[0] = (3.0 / (x_[1] - x_[0]))
               * ((y_[1] - y_[0]) / (x_[1] - x_[0]) - yp1);
    }

    for (int i = 1; i < n - 1; i++)
    {
        assert(x_[i + 1] > x_[i]);
        double sig  = (x_[i] - x_[i - 1]) / (x_[i + 1] - x_[i - 1]);
        double p    = sig * y2_[i - 1] + 2.0;
        double invp = 1. / p;
        y2_[i]      = (sig - 1.0) * invp;
        u[i]        = (6.0
                       * ((y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i])
                           - (y_[i] - y_[i - 1]) / (x_[i] - x_[i - 1]))
                       / (x_[i + 1] - x_[i - 1])
                   - sig * u[i - 1])
               * invp;
    }

    double qn, un;
    if (ypn >= 1.e30)
    {
        qn = 0.0;
        un = 0.0;
    }
    else
    {
        qn = 0.5;
        un = (3.0 / (x_[n - 1] - x_[n - 2]))
             * (ypn - (y_[n - 1] - y_[n - 2]) / (x_[n - 1] - x_[n - 2]));
    }

    y2_[n - 1] = (un - qn * u[n - 2]) / (qn * y2_[n - 2] + 1.0);

    for (int k = n - 2; k >= 0; k--)
    {
        y2_[k] = y2_[k] * y2_[k + 1] + u[k];
    }
}

double RadialSpline::splint(const double x)
{
    const int n = (int)x_.size();

    // bisection high and low
    int klo = 0;
    int khi = n - 1;

    while (khi - klo > 1)
    {
        int k = (khi + klo) / 2;
        if (x_[k] > x)
            khi = k;
        else
            klo = k;
    }

    double h = x_[khi] - x_[klo];
    assert(h > 0.);

    double a = (x_[khi] - x) / h;
    double b = (x - x_[klo]) / h;

    double y
        = a * y_[klo] + b * y_[khi]
          + h * h * (1.0 / 6.0)
                * ((a * a * a - a) * y2_[klo] + (b * b * b - b) * y2_[khi]);
    return y;
}
