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
// LDA Exchange-correlation energy and potential
// Ceperley & Alder, parametrized by Perdew and Zunger
//
////////////////////////////////////////////////////////////////////////////////

#include "LDAFunctional.h"
#include "Control.h"
#include "MGmol_MPI.h"
#include "mputils.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

const double onethird   = 1. / 3.;
const double twothird   = 2. / 3.;
const double fourthird  = 4. / 3.;
const double sevensixth = 7. / 6.;

#ifdef CRAY_T3E
double cbrt(double alpha) { return pow(alpha, 0.3333333333333333333333333); }
#endif

void LDAFunctional::computeXC(void)
{
    if (nspin_ == 1)
    {
        assert(prho_ != nullptr);
        assert(pexc_ != nullptr);
        assert(pvxc1_ != nullptr);
        for (int ir = 0; ir < np_; ir++)
        {
            xc_unpolarized(prho_[ir], pexc_[ir], pvxc1_[ir]);
        }
    }
    else
    {
        // spin polarized
        assert(prho_up_ != nullptr);
        assert(prho_dn_ != nullptr);
        assert(pexc_ != nullptr);
        assert(pvxc1_up_ != nullptr);
        assert(pvxc1_dn_ != nullptr);
        const double fz_prefac  = 1.0 / (cbrt(2.0) * 2.0 - 2.0);
        const double dfz_prefac = (fourthird)*fz_prefac;
        for (int ir = 0; ir < np_; ir++)
        {
            double excir = 0.0;
            double v_up  = 0.0;
            double v_dn  = 0.0;

            double roe_up = prho_up_[ir];
            double roe_dn = prho_dn_[ir];
            double roe    = roe_up + roe_dn;

            if (roe > 0.0)
            {
                double zeta = (roe_up - roe_dn) / roe;

                double zp1    = 1.0 + zeta;
                double zm1    = 1.0 - zeta;
                double zp1_13 = cbrt(zp1);
                double zm1_13 = cbrt(zm1);
                double fz     = fz_prefac * (zp1_13 * zp1 + zm1_13 * zm1 - 2.0);
                double dfz    = dfz_prefac * (zp1_13 - zm1_13);

                POTDTYPE xc_u, xc_p, v_u, v_p;
                xc_unpolarized(roe, xc_u, v_u);
                xc_polarized(roe, xc_p, v_p);

                double xc_pu = (double)xc_p - (double)xc_u;
                excir        = (double)xc_u + fz * xc_pu;

                double v = (double)v_u + fz * ((double)v_p - (double)v_u);
                v_up     = v + xc_pu * (1.0 - zeta) * dfz;
                v_dn     = v + xc_pu * (-1.0 - zeta) * dfz;
            }

            pvxc1_up_[ir] = (POTDTYPE)v_up;
            pvxc1_dn_[ir] = (POTDTYPE)v_dn;
            pexc_[ir]     = (POTDTYPE)excir;
        }
    }
}

double LDAFunctional::computeRhoDotExc() const
{
    double exc = 0.;
    if (nspin_ == 1)
    {
        exc = LinearAlgebraUtils<MemorySpace::Host>::MPdot(
            np_, prho_, &exc_[0]);
    }
    else
    {
        double exc_temp = 0.;
        for (int i = 0; i < np_; i++)
        {
            const double rh = prho_up_[i] + prho_dn_[i];
            exc_temp += rh * (double)pexc_[i];
        }
        exc += (POTDTYPE)exc_temp;
    }
    double sum      = 0.;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int rc          = mmpi.allreduce(&exc, &sum, 1, MPI_SUM);
    if (rc != MPI_SUCCESS)
    {
        (*MPIdata::sout) << "MPI_Allreduce double sum failed!!!" << endl;
        mmpi.abort();
    }
    exc = (POTDTYPE)sum;
    return exc;
}

void LDAFunctional::xc_unpolarized(
    const RHODTYPE rh, POTDTYPE& ee, POTDTYPE& vv)
{
    // compute LDA xc energy and potential, unpolarized
    // const double third=1.0/3.0;
    // c1 is (3.D0/(4.D0*pi))**third
    const double c1 = 0.6203504908994001;
    // alpha = (4/(9*pi))**third = 0.521061761198
    // const double alpha = 0.521061761198;
    // c2 = -(3/(4*pi)) / alpha = -0.458165293283
    // const double c2 = -0.458165293283;
    // c3 = (4/3) * c2 = -0.610887057711
    const double c3 = -0.610887057711;

    const double A  = 0.0311;
    const double B  = -0.048;
    const double b1 = 1.0529;
    const double b2 = 0.3334;
    const double G  = -0.1423;

    // C from the PZ paper: const double C  =  0.0020;
    // D from the PZ paper: const double D  = -0.0116;
    // C and D by matching Ec and Vc at rs=1
    const double D = G / (1.0 + b1 + b2) - B;
    const double C
        = -A - D - G * ((b1 / 2.0 + b2) / ((1.0 + b1 + b2) * (1.0 + b1 + b2)));

    if (rh > 0.0)
    {
        const double ro13 = cbrt(rh);
        const double rs   = c1 / ro13;

        // Next line : exchange part in Hartree units
        const double vx = c3 / rs;
        const double ex = 0.75 * vx;

        // Next lines : Ceperley & Alder correlation (Zunger & Perdew)
        double ec = 0.0, vc = 0.0;
        if (rs < 1.0)
        {
            const double logrs = log(rs);
            ec                 = A * logrs + B + C * rs * logrs + D * rs;
            vc = A * logrs + (B - A * onethird) + (twothird)*C * rs * logrs
                 + ((2.0 * D - C) * onethird) * rs;
        }
        else
        {
            const double sqrtrs  = sqrt(rs);
            const double den     = 1.0 + b1 * sqrtrs + b2 * rs;
            const double inv_den = 1. / den;
            ec                   = G * inv_den;
            vc = ec * (1.0 + sevensixth * b1 * sqrtrs + fourthird * b2 * rs)
                 * inv_den;
        }
        ee = (POTDTYPE)(ex + ec);
        vv = (POTDTYPE)(vx + vc);
    }
    else
    {
        ee = 0.;
        vv = 0.;
    }
}

void LDAFunctional::xc_polarized(const RHODTYPE rh, POTDTYPE& ee, POTDTYPE& vv)
{
    // compute LDA polarized XC energy and potential

    // const double third=1.0/3.0;
    // c1 is (3.D0/(4.D0*pi))**third
    const double c1 = 0.6203504908994001;
    // alpha = (4/(9*pi))**third = 0.521061761198
    // const double alpha = 0.521061761198;
    // c2 = -(3/(4*pi)) / alpha = -0.458165293283
    // const double c2 = -0.458165293283;
    // c3 = (4/3) * c2 = -0.610887057711
    // const double c3 = -0.610887057711;
    // c4 = 2**third * c3
    const double c4 = -0.769669463118;

    const double A  = 0.01555;
    const double B  = -0.0269;
    const double b1 = 1.3981;
    const double b2 = 0.2611;
    const double G  = -0.0843;
    // C from PZ paper: const double C   =  0.0007;
    // D from PZ paper: const double D   = -0.0048;
    // C and D by matching Ec and Vc at rs=1
    const double D = G / (1.0 + b1 + b2) - B;
    const double C
        = -A - D - G * ((b1 / 2.0 + b2) / ((1.0 + b1 + b2) * (1.0 + b1 + b2)));

    ee = 0.0;
    vv = 0.0;

    if (rh > 0.0)
    {
        double ro13 = cbrt(rh);
        double rs   = c1 / ro13;

        double ex = 0.0, vx = 0.0, ec = 0.0, vc = 0.0;

        // Next line : exchange part in Hartree units
        vx = c4 / rs;
        ex = 0.75 * vx;

        // Next lines : Ceperley & Alder correlation (Zunger & Perdew)
        if (rs < 1.0)
        {
            double logrs = log(rs);
            ec           = A * logrs + B + C * rs * logrs + D * rs;
            vc = A * logrs + (B - A * onethird) + (twothird)*C * rs * logrs
                 + ((2.0 * D - C) * onethird) * rs;
        }
        else
        {
            double sqrtrs = sqrt(rs);
            double den    = 1.0 + b1 * sqrtrs + b2 * rs;
            ec            = G / den;
            vc = ec * (1.0 + sevensixth * b1 * sqrtrs + (fourthird)*b2 * rs)
                 / den;
        }
        ee = (POTDTYPE)(ex + ec);
        vv = (POTDTYPE)(vx + vc);
    }
}
