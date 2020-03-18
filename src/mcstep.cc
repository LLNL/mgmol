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
// mcstep.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id$

#include "mcstep.h"
#include <cmath>
#include <vector>

int mcstep(double& stx, double& fx, double& dx, double& sty, double& fy,
    double& dy, double& stp, const double fp, const double dp, bool& brackt,
    const double stpmin, const double stpmax)
{
    // Adapted from fortran77 subroutine mcstep
    // by Jorge J. More', David J. Thuente
    // Argonne National Laboratory. Minpack project. June 1983
    //
    // the purpose of mcstep is to compute a safeguarded step for
    // a linesearch and to update an interval of uncertainty for
    // a minimizer of the function.
    //
    // the parameter stx contains the step with the least function
    // value. the parameter stp contains the current step. it is
    // assumed that the derivative at stx is negative in the
    // direction of the step. if brackt is set true then a
    // minimizer has been bracketed in an interval of uncertainty
    // with endpoints stx and sty.
    //
    //   stx, fx, and dx are variables which specify the step,
    //     the function, and the derivative at the best step obtained
    //     so far. the derivative must be negative in the direction
    //     of the step, that is, dx and stp-stx must have opposite
    //     signs. on output these parameters are updated appropriately.
    //
    //   sty, fy, and dy are variables which specify the step,
    //     the function, and the derivative at the other endpoint of
    //     the interval of uncertainty. on output these parameters are
    //     updated appropriately.
    //
    //   stp, fp, and dp are variables which specify the step,
    //     the function, and the derivative at the current step.
    //     if brackt is set true then on input stp must be
    //     between stx and sty. on output stp is set to the new step.
    //
    //   brackt is a boolean which specifies if a minimizer
    //     has been bracketed. if the minimizer has not been bracketed
    //     then on input brackt must be set false. if the minimizer
    //     is bracketed then on output brackt is set true.
    //
    //   stpmin and stpmax are input variables which specify lower
    //     and upper bounds for the step.
    //
    //   mcstep returns info, an integer variable set as follows:
    //     if info = 1,2,3,4,5, then the step has been computed
    //     according to one of the five cases below. otherwise
    //     info = 0, and this indicates improper input parameters.

    double stpf;
    int info = 0;
    bool bound;

    // check the input parameters for errors.
    if ((brackt && (stp <= std::min(stx, sty) || stp >= std::max(stx, sty)))
        || dx * (stp - stx) >= 0. || stpmax < stpmin)
        return info;

    // determine if the derivatives have opposite sign.
    double sgnd = dp * (dx / std::abs(dx));

    // first case. a higher function value.
    // the minimum is bracketed. if the cubic step is closer
    // to stx than the quadratic step, the cubic step is taken,
    // else the average of the cubic and quadratic steps is taken.
    if (fp > fx)
    {
        info         = 1;
        bound        = true;
        double theta = 3. * (fx - fp) / (stp - stx) + dx + dp;
        double s     = std::max(std::abs(theta), std::abs(dx));
        s            = std::max(s, std::abs(dp));
        double invs  = 1. / s;
        double ts    = theta * invs;
        double gamma = s * sqrt(ts * ts - (dx * invs) * (dp * invs));
        if (stp < stx) gamma = -gamma;
        double p    = (gamma - dx) + theta;
        double q    = ((gamma - dx) + gamma) + dp;
        double r    = p / q;
        double stpc = stx + r * (stp - stx);
        double stpq
            = stx + (0.5 * (dx / ((fx - fp) / (stp - stx) + dx))) * (stp - stx);
        if (std::abs(stpc - stx) < std::abs(stpq - stx))
        {
            stpf = stpc;
        }
        else
        {
            stpf = stpc + 0.5 * (stpq - stpc);
        }
        brackt = true;

        // second case. a lower function value and derivatives of
        // opposite sign. the minimum is bracketed. if the cubic
        // step is closer to stx than the quadratic (secant) step,
        // the cubic step is taken, else the quadratic step is taken.
    }
    else if (sgnd < 0.0)
    {
        info         = 2;
        bound        = false;
        double theta = 3. * (fx - fp) / (stp - stx) + dx + dp;
        double s     = std::max(std::abs(theta), std::abs(dx));
        s            = std::max(s, std::abs(dp));
        double invs  = 1. / s;
        double ts    = theta * invs;
        double gamma = s * sqrt(ts * ts - (dx * invs) * (dp * invs));
        if (stp > stx) gamma = -gamma;
        double p    = (gamma - dp) + theta;
        double q    = ((gamma - dp) + gamma) + dx;
        double r    = p / q;
        double stpc = stp + r * (stx - stp);
        double stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if (std::abs(stpc - stp) > std::abs(stpq - stp))
        {
            stpf = stpc;
        }
        else
        {
            stpf = stpq;
        }
        brackt = true;

        // third case. a lower function value, derivatives of the
        // same sign, and the magnitude of the derivative decreases.
        // the cubic step is only used if the cubic tends to infinity
        // in the direction of the step or if the minimum of the cubic
        // is beyond stp. otherwise the cubic step is defined to be
        // either stpmin or stpmax. the quadratic (secant) step is also
        // computed and if the minimum is bracketed then the the step
        // closest to stx is taken, else the step farthest away is taken.
    }
    else if (std::abs(dp) < std::abs(dx))
    {
        info         = 3;
        bound        = true;
        double theta = 3. * (fx - fp) / (stp - stx) + dx + dp;
        double s     = std::max(std::abs(theta), std::abs(dx));
        s            = std::max(s, std::abs(dp));

        // the case gamma = 0 only arises if the cubic does not tend
        // to infinity in the direction of the step.
        double invs = 1. / s;
        double ts   = theta * invs;
        double gamma
            = s * std::sqrt(std::max(0., ts * ts - (dx * invs) * (dp * invs)));
        if (stp > stx) gamma = -gamma;
        double p = (gamma - dp) + theta;
        double q = (gamma + (dx - dp)) + gamma;
        double r = p / q;
        double stpc;
        if (r < 0. && gamma != 0.)
        {
            stpc = stp + r * (stx - stp);
        }
        else if (stp > stx)
        {
            stpc = stpmax;
        }
        else
        {
            stpc = stpmin;
        }
        double stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if (brackt)
        {
            if (std::abs(stp - stpc) < std::abs(stp - stpq))
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpq;
            }
        }
        else
        {
            if (std::abs(stp - stpc) > std::abs(stp - stpq))
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpq;
            }
        }
        // fourth case. a lower function value, derivatives of the
        // same sign, and the magnitude of the derivative does
        // not decrease. if the minimum is not bracketed, the step
        // is either stpmin or stpmax, else the cubic step is taken.
    }
    else
    {
        info  = 4;
        bound = false;
        if (brackt)
        {
            double theta = 3. * (fp - fy) / (sty - stp) + dy + dp;
            double s     = std::max(std::abs(theta), std::abs(dy));
            s            = std::max(s, std::abs(dp));
            double invs  = 1. / s;
            double ts    = theta * invs;
            double gamma = s * sqrt(ts * ts - (dy * invs) * (dp * invs));
            if (stp > sty) gamma = -gamma;
            double p    = (gamma - dp) + theta;
            double q    = ((gamma - dp) + gamma) + dy;
            double r    = p / q;
            double stpc = stp + r * (sty - stp);
            stpf        = stpc;
        }
        else if (stp > stx)
        {
            stpf = stpmax;
        }
        else
        {
            stpf = stpmin;
        }
    }

    // update the interval of uncertainty. this update does not
    // depend on the new step or the case analysis above.
    if (fp > fx)
    {
        sty = stp;
        fy  = fp;
        dy  = dp;
    }
    else
    {
        if (sgnd < 0.)
        {
            sty = stx;
            fy  = fx;
            dy  = dx;
        }
        stx = stp;
        fx  = fp;
        dx  = dp;
    }

    // compute the new step and safeguard it.
    stpf = std::min(stpmax, stpf);
    stpf = std::max(stpmin, stpf);
    stp  = stpf;
    if (brackt && bound)
    {
        if (sty > stx)
        {
            stp = std::min(stx + 0.66 * (sty - stx), stp);
        }
        else
        {
            stp = std::max(stx + 0.66 * (sty - stx), stp);
        }
    }
    return info;
}
