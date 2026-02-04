// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_CHEBYSHEV_APPROX_INTERFACE_H
#define MGMOL_CHEBYSHEV_APPROX_INTERFACE_H

#include "ChebyshevApproximationFunction.h"
#include "LocalMatrices.h"

#include <vector>

class ChebyshevApproximationInterface
{
private:
    static Timer compute_coeffs_tm_;
    static Timer double_loop_tm_;

protected:
    int max_order_;
    int order_; // order of approximation. Also size of interpolation points for
    double extents_[2]; // interval where function is defined
    std::vector<double> interp_points_; // Chebyshev interpolation points

    LocalMatrices<double, MemorySpace::Host>*
        cmat_; // matrix to hold cosine information for
               // generating coefficients coeffs
    ChebyshevApproximationFunction*
        chebfunc_; // Pointer to function to be approximated
    std::vector<double>
        angles_; // store "angle" data for computing Chebyshev nodes
    std::vector<double> coeffs_; // store Chebyshev coefficients

    // compute Chebyshev interpolation points mapped to the interval [a, b]
    void computeInterpolationPoints();

    // evaluate function at interpolation points
    void evaluateFunction(
        const std::vector<double>& fnodes, std::vector<double>& fvals);

    // shift and scale arbitrary nodes to Chebyshev interval [-1, 1]
    void scalePointsToChebyshevInterval(
        const std::vector<double>& points, std::vector<double>& scaled_points);

public:
    // constructor
    ChebyshevApproximationInterface() {};

    // build the Chebyshev coefficients for the interval [a, b]
    void computeChebyshevCoeffs();

    // Compute Chebyshev approximation given function points.
    // Assume that points are sorted from left to right on the real axis.
    // Example: (sorted) eigenvalues of some matrix.
    void computeChebyshevApproximation(
        std::vector<double>& points, std::vector<double>& vals);

    int order() { return order_; }
    void resetOrder(const int order)
    {
        max_order_ = order;
        order_     = order;
    }
    virtual ~ChebyshevApproximationInterface() {};

    static void printTimers(std::ostream& os)
    {
        compute_coeffs_tm_.print(os);
        double_loop_tm_.print(os);
    }
};

#endif
