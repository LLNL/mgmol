// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ChebyshevApproximationInterface.h"

Timer ChebyshevApproximationInterface::compute_coeffs_tm_(
    "ChebyshevApproximationInterface::compChebCoeffs");
Timer ChebyshevApproximationInterface::double_loop_tm_(
    "ChebyshevApproximationInterface::double_loop");

void ChebyshevApproximationInterface::evaluateFunction(
    const std::vector<double>& fnodes, std::vector<double>& fvals)
{
    fvals.clear();
    fvals = chebfunc_->eval(fnodes);
}

void ChebyshevApproximationInterface::computeChebyshevCoeffs()
{
    compute_coeffs_tm_.start();
    order_ = max_order_;
    int n  = order_;
    // define scaling variable here
    const double fac = 2.0 / static_cast<double>(n);
    // evaluate function at interpolation points
    std::vector<double> fvals;
    evaluateFunction(interp_points_, fvals);

    // build coefficients
    coeffs_.clear();
    coeffs_.reserve(order_);

    // initialize first entry of coeffs_
    double fsum = 0.;
    for (const double& value : fvals)
        fsum += value;
    // scale and store first coeff
    fsum /= static_cast<double>(n);
    coeffs_.push_back(fsum);

    // compute remaining coeffs_ values
    double_loop_tm_.start();

    LocalVector<double, MemorySpace::Host> fvec(fvals);
    coeffs_.resize(n, 0.);
    LocalVector<double, MemorySpace::Host> cvec(coeffs_);

    cmat_->matvec(fvec, cvec);
    (cvec.data())[0] = fsum;
    cvec.swap(coeffs_);

    double_loop_tm_.stop();

    // scale only last n-1 coefficients
    const int n1 = n - 1;
    LinearAlgebraUtils<MemorySpace::Host>::MPscal(n1, fac, &coeffs_[1]);

    compute_coeffs_tm_.stop();
}

// Compute Chebyshev approximation given function points.
// Assume that points are sorted from left to right on the real axis.
// Example: (sorted) eigenvalues of some matrix.
void ChebyshevApproximationInterface::computeChebyshevApproximation(
    std::vector<double>& points, std::vector<double>& vals)
{
    assert(static_cast<int>(points.size()) == order_);
    vals.clear();
    // scale points to be in the range [-1, 1]
    std::vector<double> scaled_points;
    scalePointsToChebyshevInterval(points, scaled_points);

    for (int j = 0; j < order_; j++)
        assert(std::abs(scaled_points[j]) <= 1.);

    // compute Chebyshev coefficients
    computeChebyshevCoeffs();

    // initialize vals
    for (int j = 0; j < order_; j++)
        vals.push_back(coeffs_[0]);

    // loop to compute approximation values
    for (int k = 1; k < order_; k++)
    {
        for (int j = 0; j < order_; j++)
        {
            double prod = k * std::acos(scaled_points[j]);
            double t_kj = std::cos(prod);
            vals[j] += coeffs_[k] * t_kj;
        }
    }

    for (int j = 0; j < order_; j++)
        assert(vals[j] == vals[j]);
}

void ChebyshevApproximationInterface::computeInterpolationPoints()
{
    assert(order_ > 0);

    interp_points_.clear();
    angles_.clear();
    const double a  = extents_[0];
    const double b  = extents_[1];
    const double a1 = -1.;
    const double b1 = 1.;

    const double alpha = M_PI / (2. * order_);
    const double beta  = (b - a) / (b1 - a1);
    for (int i = 0; i < order_; i++)
    {
        double ang = (2. * (i + 1.) - 1.) * alpha;
        angles_.push_back(ang);
        double x_i = std::cos(ang);
        // map node to interval [a, b]
        const double scaled_node = a + beta * (x_i - a1);

        assert(scaled_node == scaled_node);
        interp_points_.push_back(scaled_node);
    }

    double_loop_tm_.start();
    if (cmat_) delete cmat_;
    cmat_ = new LocalMatrices<double, MemorySpace::Host>(1, order_, order_);
    std::vector<double> tmp(order_ * order_);
    for (int i = 1; i < order_; i++)
    {
        for (int k = 0; k < order_; k++)
        {
            double iang         = i * angles_[k];
            double val          = std::cos(iang);
            tmp[i + k * order_] = val;
        }
    }
    cmat_->setValues(tmp.data(), order_);

    double_loop_tm_.stop();
}

void ChebyshevApproximationInterface::scalePointsToChebyshevInterval(
    const std::vector<double>& points, std::vector<double>& scaled_points)
{
    const double a1 = -1;
    const double b1 = 1;
    // get extents of points
    const double a = points[0];
    const double b = points[order_ - 1];

    const double alpha = (b1 - a1) / (b - a);

    scaled_points.clear();
    scaled_points.reserve(points.size());

    for (const double& val : points)
    {
        double sp = a1 + alpha * (val - a);
        assert(sp == sp);

        scaled_points.push_back(sp);
    }
}
