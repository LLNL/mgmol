// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cassert>
#include <iostream>
#include <math.h>

#include "radial_functions.h"

const double e0 = 78.36;
// const double  e0 = 80.;
// const double  e0 = 1.;

#ifndef NDEBUG
const double tol = 1.e-15;
#endif

namespace pb
{

// nu2, nu4, nu6 from L. David and M.J. Field, JCC 18 (1997)
double nu2(const double r, const double rcav, const double rc)
{
    double out;

    const double s  = (r - rcav) / rc;
    const double s2 = s * s;

    if (s > 0.)
    {
        out = 1. - (1. + s2) * exp(-s2);
    }
    else
    {
        out = 0.;
    }
    return out;
}

double nu4(const double r, const double rcav, const double rc)
{
    double out;
    const double s  = (r - rcav) / rc;
    const double s2 = s * s;
    const double s4 = s2 * s2;
    if (s > 0.)
    {
        out = 1. - (1. + s2 + 0.5 * s4) * exp(-s2);
    }
    else
    {
        out = 0.;
    }

    assert(out + tol >= 0.);
    return std::max(out, 0.);
}

double nu6(const double r, const double rcav, const double rc)
{
    double out;

    assert(rc > 0.);

    const double s  = (r - rcav) / rc;
    const double s2 = s * s;
    const double s4 = s2 * s2;
    const double s6 = s2 * s4;

    if (s > 0.)
    {
        out = 1. - (1. + s2 + 0.5 * s4 + 0.1666666666666 * s6) * exp(-s2);
        // out = 1.-(1.+s2+0.5*s4+0.2*s6)*exp(-s2);
    }
    else
    {
        out = 0.;
    }
    return out;
}

double epsilon2(const double r, const double a, const double rc)
{
    double out = 1. + (e0 - 1.) * nu2(r, a, rc);

    return out;
}

double epsilon4(const double r, const double a, const double rc)
{
    assert(r >= 0.);

    double out = 1. + (e0 - 1.) * nu4(r, a, rc);

    return out;
}

double epsilon4(const double r, const double a)
{
    assert(r >= 0.);

    double out = 1. + (e0 - 1.) * nu4(r, a, 1.);

    assert(out + tol >= 1.);
    return out;
}

double epsilon6(const double r, const double a, const double rc)
{
    double out = 1. + (e0 - 1.) * nu4(r, a, rc);

    return out;
}

// normalized Gaussian
double gaussian(const double r, const double rc)
{
    assert(rc > 0.);

    return exp(-(r * r) / (2. * rc * rc))
           / (rc * rc * rc * 2. * M_PI * sqrt(2. * M_PI));
}

double comp_potential(const double r, const double rc)
{
    assert(rc > 0.);

    if (r > 0.00001)
        return (erf(2. * r / rc) - erf(r / rc)) / (M_PI * 4. * r);
    else
        return 1. / (2. * rc * sqrt(M_PI) * M_PI);
}

// solution of -Lap u=f for f Gaussian
double potential_of_gaussian(const double r, const double rc)
{
    assert(rc > 0.);

    if (r > 0.00001)
        return erf(r / rc) / r;
    else
        return 2. / (rc * sqrt(M_PI));
}

double comp_charge(const double r, const double rc)
{
    assert(rc > 0.);

    const double r2 = r * r / (rc * rc);

    return (8. * exp(-4. * r2) - exp(-r2)) / (rc * rc * rc * M_PI * sqrt(M_PI));
}

double rho2(const double r, const double a)
{
    double alpha;

    const double r2 = r * r;

    if (r > a)
    {
        const double ar2 = (a - r) * (a - r);

        alpha
            = (4 * (-8 + exp(3 * r2))
                  * (1
                      + (-1 + e0)
                            * (1
                                - (1 + pow(a, 2) - 2 * a * r + r2) / exp(ar2))))
                  / exp(4 * r2)
              - ((-1 + e0)
                    * ((2 * (a - r)) / exp(ar2)
                        + (2 * (1 + ar2 + 0. * pow(a - r, 4)) * (-a + r))
                              / exp(ar2))
                    * ((2 * (-2 + exp(3 * r2)) * r) / exp(4 * r2)
                        + sqrt(M_PI) * (erf(2. * r) - erf(r))))
                    / r2;
    }
    else
    {
        alpha = 4. * (-8. * exp(-4. * r2) + exp(-r2));
    }

    return -alpha;
}
double rho4(const double r, const double a)
{
    double alpha;

    const double r2 = r * r;

    if (r > (a + 9)) return 0.;

    if (r > a)
    {

        const double ar2 = (a - r) * (a - r);
        const double ar3 = ar2 * (a - r);
        const double ar4 = ar3 * (a - r);
        const double ar5 = ar4 * (a - r);

        alpha = ((4 * (-8 + exp(3 * r2))
                     * (1 + (-1 + e0) * (1 - (1 + ar2 + 0.5 * ar4) / exp(ar2)))
                     * r2)
                        / exp(4 * r2)
                    + exp(-a * a + 2 * a * r - 5 * r2) * (-1 + e0) * ar5
                          * (2 * (-2 + exp(3 * r2)) * r
                              + exp(4 * r2) * sqrt(M_PI)
                                    * (erf(2. * r) - erf(r))))
                / (4. * pow(M_PI, 1.5) * r2);
    }
    else
    {
        alpha = (-8 + exp(3 * r2)) / (exp(4 * r2) * sqrt(M_PI) * M_PI);
    }

    return -alpha;
}

} // namespace pb
