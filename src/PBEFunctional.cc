// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "PBEFunctional.h"
#include "Control.h"
#include "MGmol_MPI.h"
#include "mputils.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

const static double uk = 0.804;

PBEFunctional::PBEFunctional(std::vector<std::vector<RHODTYPE>>& rhoe)
    : XCFunctional(rhoe)
{
    pgrad_rho_[0] = pgrad_rho_[1] = pgrad_rho_[2] = nullptr;
    pgrad_rho_up_[0] = pgrad_rho_up_[1] = pgrad_rho_up_[2] = nullptr;
    pgrad_rho_dn_[0] = pgrad_rho_dn_[1] = pgrad_rho_dn_[2] = nullptr;
    if (nspin_ == 1)
    {
        exc_.resize(np_);
        vxc1_.resize(np_);
        vxc2_.resize(np_);
        grad_rho_[0].resize(np_);
        grad_rho_[1].resize(np_);
        grad_rho_[2].resize(np_);
        pgrad_rho_[0] = &grad_rho_[0][0];
        pgrad_rho_[1] = &grad_rho_[1][0];
        pgrad_rho_[2] = &grad_rho_[2][0];
        pexc_         = &exc_[0];
        pvxc1_        = &vxc1_[0];
        pvxc2_        = &vxc2_[0];
    }
    else
    {
        exc_up_.resize(np_);
        exc_dn_.resize(np_);
        vxc1_up_.resize(np_);
        vxc1_dn_.resize(np_);
        vxc2_upup_.resize(np_);
        vxc2_updn_.resize(np_);
        vxc2_dnup_.resize(np_);
        vxc2_dndn_.resize(np_);
        grad_rho_up_[0].resize(np_);
        grad_rho_up_[1].resize(np_);
        grad_rho_up_[2].resize(np_);
        grad_rho_dn_[0].resize(np_);
        grad_rho_dn_[1].resize(np_);
        grad_rho_dn_[2].resize(np_);

        pgrad_rho_up_[0] = &grad_rho_up_[0][0];
        pgrad_rho_up_[1] = &grad_rho_up_[1][0];
        pgrad_rho_up_[2] = &grad_rho_up_[2][0];
        pgrad_rho_dn_[0] = &grad_rho_dn_[0][0];
        pgrad_rho_dn_[1] = &grad_rho_dn_[1][0];
        pgrad_rho_dn_[2] = &grad_rho_dn_[2][0];
        pexc_up_         = &exc_up_[0];
        pexc_dn_         = &exc_dn_[0];
        pvxc1_up_        = &vxc1_up_[0];
        pvxc1_dn_        = &vxc1_dn_[0];
        pvxc2_upup_      = &vxc2_upup_[0];
        pvxc2_updn_      = &vxc2_updn_[0];
        pvxc2_dnup_      = &vxc2_dnup_[0];
        pvxc2_dndn_      = &vxc2_dndn_[0];
    }
}

void PBEFunctional::computeXC(void)
{
    if (nspin_ == 1)
    {
        assert(prho_ != nullptr);
        assert(pgrad_rho_[0] != nullptr && pgrad_rho_[1] != nullptr
               && pgrad_rho_[2] != nullptr);
        assert(pexc_ != nullptr);
        assert(pvxc1_ != nullptr);
        assert(pvxc2_ != nullptr);

        const RHODTYPE* const grad0 = pgrad_rho_[0];
        const RHODTYPE* const grad1 = pgrad_rho_[1];
        const RHODTYPE* const grad2 = pgrad_rho_[2];
        for (int i = 0; i < np_; i++)
        {
            double gp = (double)grad0[i] * (double)grad0[i]
                        + (double)grad1[i] * (double)grad1[i]
                        + (double)grad2[i] * (double)grad2[i];
            RHODTYPE grad = (POTDTYPE)sqrt(gp);

            // compute energy density pexc_ and potentials components
            // pvxc1_ and pvxc2_
            excpbe(prho_[i], grad, pexc_ + i, pvxc1_ + i, pvxc2_ + i);
        }
    }
    else
    {
        assert(prho_up_ != nullptr);
        assert(prho_dn_ != nullptr);
        assert(pgrad_rho_up_[0] != nullptr && pgrad_rho_up_[1] != nullptr
               && pgrad_rho_up_[2] != nullptr);
        assert(pgrad_rho_dn_[0] != nullptr && pgrad_rho_dn_[1] != nullptr
               && pgrad_rho_dn_[2] != nullptr);
        assert(pexc_up_ != nullptr);
        assert(pexc_dn_ != nullptr);
        assert(pvxc1_up_ != nullptr);
        assert(pvxc1_dn_ != nullptr);
        assert(pvxc2_upup_ != nullptr);
        assert(pvxc2_updn_ != nullptr);
        assert(pvxc2_dnup_ != nullptr);
        assert(pvxc2_dndn_ != nullptr);

        for (int i = 0; i < np_; i++)
        {
            RHODTYPE grx_up = pgrad_rho_up_[0][i];
            RHODTYPE gry_up = pgrad_rho_up_[1][i];
            RHODTYPE grz_up = pgrad_rho_up_[2][i];
            RHODTYPE grx_dn = pgrad_rho_dn_[0][i];
            RHODTYPE gry_dn = pgrad_rho_dn_[1][i];
            RHODTYPE grz_dn = pgrad_rho_dn_[2][i];
            double grx      = (double)grx_up + (double)grx_dn;
            double gry      = (double)gry_up + (double)gry_dn;
            double grz      = (double)grz_up + (double)grz_dn;
            double grxyz_up = (double)grx_up * (double)grx_up
                              + (double)gry_up * (double)gry_up
                              + (double)grz_up * (double)grz_up;
            double grad_up  = sqrt(grxyz_up);
            double grxyz_dn = (double)grx_dn * (double)grx_dn
                              + (double)gry_dn * (double)gry_dn
                              + (double)grz_dn * (double)grz_dn;
            double grad_dn = sqrt(grxyz_dn);
            double grad    = sqrt(grx * grx + gry * gry + grz * grz);
            excpbe_sp(prho_up_[i], prho_dn_[i], grad_up, grad_dn, grad,
                &pexc_up_[i], &pexc_dn_[i], &pvxc1_up_[i], &pvxc1_dn_[i],
                &pvxc2_upup_[i], &pvxc2_dndn_[i], &pvxc2_updn_[i],
                &pvxc2_dnup_[i]);
        }
    }
}

double PBEFunctional::computeRhoDotExc() const
{
    double exc = 0.;
    if (nspin_ == 1)
    {
        exc = LinearAlgebraUtils<MemorySpace::Host>::MPdot(
            np_, prho_, &exc_[0]);
    }
    else
    {
        for (int i = 0; i < np_; i++)
        {
            exc += prho_up_[i] * pexc_up_[i];
            exc += prho_dn_[i] * pexc_dn_[i];
        }
    }
    double sum      = 0.;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int rc          = mmpi.allreduce(&exc, &sum, 1, MPI_SUM);
    if (rc != MPI_SUCCESS)
    {
        (*MPIdata::sout) << "MPI_Allreduce double sum failed!!!" << std::endl;
        mmpi.abort();
    }
    exc = sum;
    return exc;
}

////////////////////////////////////////////////////////////////////////////////
//
//  excpbe: PBE exchange-correlation
//  K.Burke's modification of PW91 codes, May 14, 1996.
//  Modified again by K.Burke, June 29, 1996, with simpler Fx(s)
//  Translated into C and modified by F.Gygi, Dec 9, 1996.
//
//  input:
//    rho:  density
//    grad: abs(grad(rho))
//  output:
//    exc: exchange-correlation energy per electron
//    vxc1, vxc2 : quantities such that the total exchange potential is:
//
//      vxc = vxc1 + div ( vxc2 * grad(n) )
//
//  References:
//  [a] J.P.Perdew, K.Burke, and M.Ernzerhof,
//      "Generalized gradient approximation made simple,
//      Phys.Rev.Lett. 77, 3865, (1996).
//  [b] J.P.Perdew and Y.Wang, Phys.Rev. B33, 8800 (1986),
//      Phys.Rev. B40, 3399 (1989) (E).
//
////////////////////////////////////////////////////////////////////////////////

void PBEFunctional::excpbe(const RHODTYPE rho, const RHODTYPE grad,
    POTDTYPE* exc, POTDTYPE* vxc1, POTDTYPE* vxc2)
{
    const double third  = 1.0 / 3.0;
    const double third4 = 4.0 / 3.0;
    const double ax     = -0.7385587663820224058; /* -0.75*pow(3.0/pi,third) */
    const double um     = 0.2195149727645171;
    const double ul     = um / uk;
    const double pi32third   = 3.09366772628014; /* (3*pi^2 ) ^(1/3) */
    const double alpha       = 1.91915829267751; /* pow(9.0*pi/4.0, third)*/
    const double seven_sixth = 7.0 / 6.0;
    const double gamma       = 0.03109069086965489; /* gamma = (1-ln2)/pi^2 */
    const double bet         = 0.06672455060314922; /* see [a] (4) */
    const double delt        = bet / gamma;

    double rtrs, twoks, rs, t, h, ecrs, pon, b, b2, t2, t4, q4, q5, t6, rsthrd,
        fac, bec, q8, q9, hb, hrs, ht, vc;

    double s, s2, p0, fxpbe, fs;
    double ex, vx1, vx2, ec, vc1, vc2;

    *exc  = 0.0;
    *vxc1 = 0.0;
    *vxc2 = 0.0;

    if (rho < 1.e-18)
    {
        return;
    }

    /* exchange */

    double rh13 = cbrt(rho);

    /* LDA exchange energy density */
    double exunif = ax * rh13;

    /* Fermi wavevector  kF = ( 3 * pi^2 n )^(1/3) */
    double fk = pi32third * rh13;

    const double invfk  = 1. / fk;
    const double invrho = 1. / rho;
    s                   = 0.5 * grad * invfk * invrho;

    /* PBE enhancement factor */

    s2                 = s * s;
    p0                 = 1.0 + ul * s2;
    const double invp0 = 1. / p0;
    fxpbe              = 1.0 + uk - uk * invp0;

    ex = exunif * fxpbe;

    /* energy done, now the potential */
    /* find first derivative of Fx w.r.t the variable s. */
    /* fs = (1/s) * d Fx / d s */

    fs = 2.0 * uk * ul * invp0 * invp0;

    vx1 = third4 * exunif * (fxpbe - s2 * fs);
    vx2 = -0.25 * exunif * fs * invfk * invfk * invrho;

    /* correlation */

    /* Find LSD contributions, using [c] (10) and Table I of [c]. */
    /* ec = unpolarized LSD correlation energy */
    /* ecrs = d ec / d rs */
    /* construct ec, using [c] (8) */

    rs                    = alpha * invfk;
    twoks                 = 2.0 * M_2_SQRTPI * sqrt(fk);
    const double invtwoks = 1. / twoks;
    t                     = grad * invrho * invtwoks;

    rtrs = sqrt(rs);
    gcor2(0.0310907, 0.2137, 7.5957, 3.5876, 1.6382, 0.49294, rtrs, &ec, &ecrs);

    /* LSD potential from [c] (A1) */
    /* ecrs = d ec / d rs [c] (A2) */

    vc = ec - rs * ecrs * third;

    /* PBE correlation energy */
    /* b = A of [a] (8) */

    pon = -ec / gamma;
    b   = delt / (exp(pon) - 1.0);
    b2  = b * b;
    t2  = t * t;
    t4  = t2 * t2;
    q4  = 1.0 + b * t2;
    q5  = q4 + b2 * t4;
    h   = gamma * log(1.0 + delt * q4 * t2 / q5);

    // Energy done, now the potential, using appendix E of [b]

    t6           = t4 * t2;
    rsthrd       = rs * third;
    fac          = delt / b + 1.0;
    bec          = b2 * fac / bet;
    q8           = q5 * q5 + delt * q4 * q5 * t2;
    double invq8 = 1. / q8;
    q9           = 1.0 + 2.0 * b * t2;
    hb           = -bet * b * t6 * (2.0 + b * t2) * invq8;
    hrs          = -rsthrd * hb * bec * ecrs;
    ht           = 2.0 * bet * q9 * invq8;

    vc1 = vc + h + hrs - t2 * ht * seven_sixth;
    vc2 = -ht * invrho * invtwoks * invtwoks;

    *exc  = (POTDTYPE)(ex + ec + h);
    *vxc1 = (POTDTYPE)(vx1 + vc1);
    *vxc2 = (POTDTYPE)(vx2 + vc2);
}

////////////////////////////////////////////////////////////////////////////////

void PBEFunctional::excpbe_sp(RHODTYPE rho_up, RHODTYPE rho_dn,
    RHODTYPE grad_up, RHODTYPE grad_dn, RHODTYPE grad, POTDTYPE* exc_up,
    POTDTYPE* exc_dn, POTDTYPE* vxc1_up, POTDTYPE* vxc1_dn, POTDTYPE* vxc2_upup,
    POTDTYPE* vxc2_dndn, POTDTYPE* vxc2_updn, POTDTYPE* vxc2_dnup)
{
    const double third  = 1.0 / 3.0;
    const double third2 = 2.0 / 3.0;
    const double third4 = 4.0 / 3.0;
    const double sixthm = -1.0 / 6.0;
    const double ax     = -0.7385587663820224058; /* -0.75*pow(3.0/pi,third) */
    const double um     = 0.2195149727645171;
    const double ul     = um / uk;
    const double pi32third    = 3.09366772628014; /* (3*pi^2 ) ^(1/3) */
    const double alpha        = 1.91915829267751; /* pow(9.0*pi/4.0, third)*/
    const double seven_sixth  = 7.0 / 6.0;
    const double four_over_pi = 1.27323954473516;
    const double gam          = 0.5198420997897463; /* gam = 2^(4/3) - 2 */
    const double fzz          = 8.0 / (9.0 * gam);
    const double gamma        = 0.03109069086965489; /* gamma = (1-ln2)/pi^2 */
    const double bet          = 0.06672455060314922; /* see [a] (4) */
    const double delt         = bet / gamma;
    const double eta = 1.e-12; // small number to avoid blowup as |zeta|->1

    double eu, eurs, ep, eprs, alfm, alfrsm;
    double ex_up, ex_dn, vx1_up, vx1_dn, vx2_up, vx2_dn, ec, vc1_up, vc1_dn,
        vc2;

    *exc_up    = 0.0;
    *exc_dn    = 0.0;
    *vxc1_up   = 0.0;
    *vxc1_dn   = 0.0;
    *vxc2_upup = 0.0;
    *vxc2_updn = 0.0;
    *vxc2_dnup = 0.0;
    *vxc2_dndn = 0.0;

    if (rho_up < 1.e-18 && rho_dn < 1.e-18)
    {
        return;
    }

    /* exchange up */

    ex_up  = 0.0;
    vx1_up = 0.0;
    vx2_up = 0.0;
    if (rho_up > 1.e-18)
    {
        double tworho = 2.0 * rho_up;
        double gr     = 2.0 * grad_up;

        double rh13 = pow(tworho, third);
        /* LDA exchange energy density */
        double exunif = ax * rh13;
        /* Fermi wavevector  kF = ( 3 * pi^2 n )^(1/3) */
        double fk = pi32third * rh13;
        double s  = gr / (2.0 * fk * tworho);
        /* PBE enhancement factor */
        double s2    = s * s;
        double p0    = 1.0 + ul * s2;
        double fxpbe = 1.0 + uk - uk / p0;
        ex_up        = exunif * fxpbe;
        /* energy done, now the potential */
        /* find first derivative of Fx w.r.t the variable s. */
        /* fs = (1/s) * d Fx / d s */
        double fs = 2.0 * uk * ul / (p0 * p0);
        vx1_up    = third4 * exunif * (fxpbe - s2 * fs);
        vx2_up    = -exunif * fs / (tworho * 4.0 * fk * fk);
    }

    /* exchange dn */

    ex_dn  = 0.0;
    vx1_dn = 0.0;
    vx2_dn = 0.0;
    if (rho_dn > 1.e-18)
    {
        double tworho = 2.0 * rho_dn;
        double gr     = 2.0 * grad_dn;

        double rh13 = pow(tworho, third);
        /* LDA exchange energy density */
        double exunif = ax * rh13;
        /* Fermi wavevector  kF = ( 3 * pi^2 n )^(1/3) */
        double fk = pi32third * rh13;
        double s  = gr / (2.0 * fk * tworho);
        /* PBE enhancement factor */
        double s2    = s * s;
        double p0    = 1.0 + ul * s2;
        double fxpbe = 1.0 + uk - uk / p0;
        ex_dn        = exunif * fxpbe;
        /* energy done, now the potential */
        /* find first derivative of Fx w.r.t the variable s. */
        /* fs = (1/s) * d Fx / d s */
        double fs = 2.0 * uk * ul / (p0 * p0);
        vx1_dn    = third4 * exunif * (fxpbe - s2 * fs);
        vx2_dn    = -exunif * fs / (tworho * 4.0 * fk * fk);
    }

    /* correlation */

    // Find LSD contributions, using [c] (10) and Table I of [c].
    // eu = unpolarized LSD correlation energy
    // eurs = d eu / d rs
    // ep = fully polarized LSD correlation energy
    // eprs = d ep / d rs
    // alfm = - spin stiffness, [c] (3)
    // alfrsm = -d alpha / d rs
    // f = spin-scaling factor from [c] (9)
    // construct ec, using [c] (8)

    double rhotot = rho_up + rho_dn;

    double rh13   = pow(rhotot, third);
    double zet    = (rho_up - rho_dn) / rhotot;
    double g      = 0.5 * (pow(1.0 + zet, third2) + pow(1.0 - zet, third2));
    double fk     = pi32third * rh13;
    double rs     = alpha / fk;
    double twoksg = 2.0 * sqrt(four_over_pi * fk) * g;
    double t      = grad / (twoksg * rhotot);

    double rtrs = sqrt(rs);
    gcor2(0.0310907, 0.2137, 7.5957, 3.5876, 1.6382, 0.49294, rtrs, &eu, &eurs);
    gcor2(0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517, rtrs, &ep,
        &eprs);
    gcor2(0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671, rtrs, &alfm,
        &alfrsm);
    double z4 = zet * zet * zet * zet;
    double f  = (pow(1.0 + zet, third4) + pow(1.0 - zet, third4) - 2.0) / gam;
    ec        = eu * (1.0 - f * z4) + ep * f * z4 - alfm * f * (1.0 - z4) / fzz;

    /* LSD potential from [c] (A1) */
    /* ecrs = d ec / d rs [c] (A2) */
    double ecrs
        = eurs * (1.0 - f * z4) + eprs * f * z4 - alfrsm * f * (1.0 - z4) / fzz;
    double fz = third4 * (pow(1.0 + zet, third) - pow(1.0 - zet, third)) / gam;
    double eczet = 4.0 * (zet * zet * zet) * f * (ep - eu + alfm / fzz)
                   + fz * (z4 * ep - z4 * eu - (1.0 - z4) * alfm / fzz);
    double comm = ec - rs * ecrs * third - zet * eczet;
    vc1_up      = comm + eczet;
    vc1_dn      = comm - eczet;

    /* PBE correlation energy */
    /* b = A of [a] (8) */

    double g3  = g * g * g;
    double pon = -ec / (g3 * gamma);
    double b   = delt / (exp(pon) - 1.0);
    double b2  = b * b;
    double t2  = t * t;
    double t4  = t2 * t2;
    double q4  = 1.0 + b * t2;
    double q5  = q4 + b2 * t4;
    double h   = g3 * gamma * log(1.0 + delt * q4 * t2 / q5);

    /* Energy done, now the potential, using appendix E of [b] */

    double g4     = g3 * g;
    double t6     = t4 * t2;
    double rsthrd = rs * third;
    double gz     = (pow((1.0 + zet) * (1.0 + zet) + eta, sixthm)
                    - pow((1.0 - zet) * (1.0 - zet) + eta, sixthm))
                * third;
    double fac  = delt / b + 1.0;
    double bg   = -3.0 * b2 * ec * fac / (bet * g4);
    double bec  = b2 * fac / (bet * g3);
    double q8   = q5 * q5 + delt * q4 * q5 * t2;
    double q9   = 1.0 + 2.0 * b * t2;
    double hb   = -bet * g3 * b * t6 * (2.0 + b * t2) / q8;
    double hrs  = -rsthrd * hb * bec * ecrs;
    double hzed = 3.0 * gz * h / g + hb * (bg * gz + bec * eczet);
    double ht   = 2.0 * bet * g3 * q9 / q8;

    double ccomm = h + hrs - t2 * ht * seven_sixth;
    double pref  = hzed - gz * t2 * ht / g;

    ccomm -= pref * zet;

    vc1_up += ccomm + pref;
    vc1_dn += ccomm - pref;
    vc2 = -ht / (rhotot * twoksg * twoksg);

    *exc_up    = (POTDTYPE)(ex_up + ec + h);
    *exc_dn    = (POTDTYPE)(ex_dn + ec + h);
    *vxc1_up   = (POTDTYPE)(vx1_up + vc1_up);
    *vxc1_dn   = (POTDTYPE)(vx1_dn + vc1_dn);
    *vxc2_upup = (POTDTYPE)(2 * vx2_up + vc2);
    *vxc2_dndn = (POTDTYPE)(2 * vx2_dn + vc2);
    *vxc2_updn = (POTDTYPE)vc2;
    *vxc2_dnup = (POTDTYPE)vc2;
}

////////////////////////////////////////////////////////////////////////////////
//
//  gcor2.c: Interpolate LSD correlation energy
//  as given by (10) of Perdew & Wang, Phys Rev B45 13244 (1992)
//  Translated into C by F.Gygi, Dec 9, 1996
//
////////////////////////////////////////////////////////////////////////////////

void PBEFunctional::gcor2(double a, double a1, double b1, double b2, double b3,
    double b4, double rtrs, double* gg, double* ggrs)
{
    double q0, q1, q2, q3;
    q0    = -2.0 * a * (1.0 + a1 * rtrs * rtrs);
    q1    = 2.0 * a * rtrs * (b1 + rtrs * (b2 + rtrs * (b3 + rtrs * b4)));
    q2    = log(1.0 + 1.0 / q1);
    *gg   = q0 * q2;
    q3    = a * (b1 / rtrs + 2.0 * b2 + rtrs * (3.0 * b3 + 4.0 * b4 * rtrs));
    *ggrs = -2.0 * a * a1 * q2 - q0 * q3 / (q1 * (1.0 + q1));
}
