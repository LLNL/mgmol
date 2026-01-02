// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef CHEBYSHEV_APPROX_H
#define CHEBYSHEV_APPROX_H

#include "ChebyshevApproximationFunction.h"
#include "ChebyshevApproximationInterface.h"

#include <iostream>
#include <vector>

template <class MatrixType>
class ChebyshevApproximation : public ChebyshevApproximationInterface
{

private:
    std::vector<MatrixType> nodesTk_; // Chebyshev nodes (polynomials)

    // Shift and scale matrix H so that spectrum is in range [-1, 1].
    // NOTE: Matrix is modified on return
    void scaleMatrixToChebyshevInterval(
        const double a, const double b, MatrixType& H);

    static Timer compute_tm_;
    static Timer compute2_tm_;
    static Timer build_nodes_tm_;

public:
    ChebyshevApproximation(const double a, const double b, const int order,
        ChebyshevApproximationFunction* func); // constructor
    void setup(); // setup some data

    // build array of Chebyshev Nodes (Polynomials)
    void buildChebyshevNodes(const double a, const double b, MatrixType& H);
    // Compute Chebyshev approximation stored data (coeffs_ and nodesTk_)
    MatrixType computeChebyshevApproximation();
    // Compute Chebyshev approximation given extents of approximation and matrix
    // argument (for computing polynomials)
    MatrixType computeChebyshevApproximation(
        const MatrixType& H, const bool recompute_coeffs = true);
    // destructor
    ~ChebyshevApproximation() override
    {
        if (cmat_) delete cmat_;
    }

    static void printTimers(std::ostream& os)
    {
        compute_tm_.print(os);
        compute2_tm_.print(os);
        build_nodes_tm_.print(os);
    }
};

template <class MatrixType>
Timer ChebyshevApproximation<MatrixType>::compute_tm_(
    "ChebyshevApproximation::compute");
template <class MatrixType>
Timer ChebyshevApproximation<MatrixType>::compute2_tm_(
    "ChebyshevApproximation::compute2");
template <class MatrixType>
Timer ChebyshevApproximation<MatrixType>::build_nodes_tm_(
    "ChebyshevApproximation::build_nodes");

#endif
