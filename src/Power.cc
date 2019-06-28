// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Power.h"

#include "GramMatrix.h"
#include "mputils.h"
#include "random.h"
#include "DistVector.h"

#include <vector>

Timer Power::compute_tm_("Power::compute");
Timer Power::compute_gen_tm_("Power::compute_gen");

std::ostream* Power::os_ = &std::cout;

// compute sum of squares of elements of vector y-theta*v
double diff2(
    std::vector<double>& y, std::vector<double>& v, const MATDTYPE theta,
    const bool verbose)
{
    double diff                             = 0.;
    std::vector<double>::const_iterator it1 = v.begin();
    for (std::vector<double>::const_iterator it2 = y.begin(); it2 != y.end();
         ++it2)
    {
        double tmp = *it2 - theta * (*it1);
        diff += tmp * tmp;

        ++it1;
    }

    if ( verbose)
        std::cout << "Power method: theta=" << theta << ", diff2=" << diff << '\n';

    return diff;
}

double power(SquareLocalMatrices<double>& A, std::vector<double>& y,
    const int maxits, const double epsilon, const bool verbose)
{
    std::vector<double> v(y.size());
    double theta = 0.;

    for (int i = 0; i < maxits; i++)
    {
        double norm = Tnrm2(y.size(), &y[0]);
        Tscal(y.size(), 1. / norm, &y[0]);
        v.swap(y);
        A.matvec(v, y, 0);
        theta              = Tdot(y.size(), &y[0], &v[0]);
        const double tol = epsilon * fabs(theta);
        // do at least 2 iterations before checking for convergence
        if (i > 1)
            if (diff2(y, v, theta, verbose) <= tol * tol)
            {
                if (verbose)
                    std::cout << "Power method converge in " << i << " iterations\n";
                break;
            }
    }

    return theta;
}

// compute extreme eigenvalues of A by power method
void Power::computeEigenInterval(
    SquareLocalMatrices<double>& A, double& emin, double& emax,
    const double epsilon, const bool verbose)
{
    compute_tm_.start();

    int maxits     = 100;

    A.shift(shift_);

    double beta1 = power(A, vec1_, maxits, epsilon, verbose);
    double e1    = beta1 - shift_;

    // shift matrix and compute other extreme eigenvalue
    A.shift(-beta1);

    double beta2 = power(A, vec2_, maxits, epsilon, verbose);

    double e2 = beta2 + beta1 - shift_;

    if (e1 < e2)
    {
        emin = e1;
        emax = e2;
    }
    else
    {
        emin = e2;
        emax = e1;
    }

    // set shift to search for e1 first at next call
    shift_ = -e2;

    compute_tm_.stop();
}
