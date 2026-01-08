// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Power.h"

#include "GramMatrix.h"
#include "LocalVector.h"
#include "ReplicatedMatrix.h"
#include "ReplicatedVector.h"
#include "SquareLocalMatrices.h"
#include "mputils.h"
#include "random.h"

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrix.h"
#include "DistVector.h"
#endif

// compute sum of squares of elements of vector y-theta*v
template <class VECTOR>
double diff2(VECTOR& y, VECTOR& v, const double theta, const bool verbose)
{
    double diff = y.scaledDiff2(v, theta);

    if (verbose)
        std::cout << "Power method: theta=" << theta << ", diff2=" << diff
                  << '\n';

    return diff;
}

template <class VECTOR, class MATRIX>
double power(MATRIX& A, VECTOR& y, const int maxits, const double epsilon,
    const bool verbose)
{
    VECTOR v(A.m());
    double theta = 0.;

    for (int i = 0; i < maxits; i++)
    {
        y.normalize();
        v.swap(y);
        A.matvec(v, y);
        theta            = y.dot(v);
        const double tol = epsilon * fabs(theta);
        // do at least 2 iterations before checking for convergence
        if (i > 1)
            if (diff2(y, v, theta, verbose) <= tol * tol)
            {
                if (verbose)
                    std::cout << "Power method converge in " << i
                              << " iterations\n";
                break;
            }
    }

    return theta;
}

// compute extreme eigenvalues of A by power method
template <class VECTOR, class MATRIX>
void Power<VECTOR, MATRIX>::computeEigenInterval(const MATRIX& A, double& emin,
    double& emax, const double epsilon, const bool verbose)
{
    compute_tm_.start();

    int maxits = 100;

    MATRIX B(A);
    B.shift(shift_);

    double beta1 = power(B, vec1_, maxits, epsilon, verbose);
    double e1    = beta1 - shift_;

    // shift matrix and compute other extreme eigenvalue
    B.shift(-beta1);

    double beta2 = power(B, vec2_, maxits, epsilon, verbose);

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

template class Power<LocalVector<double, MemorySpace::Host>,
    SquareLocalMatrices<double, MemorySpace::Host>>;
#ifdef MGMOL_USE_SCALAPACK
template class Power<dist_matrix::DistVector<double>,
    dist_matrix::DistMatrix<double>>;
#endif
// template class Power<ReplicatedVector, ReplicatedMatrix>;
