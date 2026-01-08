// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ChebyshevApproximation.h"
#include "MPIdata.h"
#include "ReplicatedMatrix.h"

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrix.h"
#endif

#include <iostream>

template <class MatrixType>
ChebyshevApproximation<MatrixType>::ChebyshevApproximation(const double a,
    const double b, const int order, ChebyshevApproximationFunction* func)
{
    assert(a < b);

    extents_[0] = a;
    extents_[1] = b;
    max_order_  = order;
    order_      = order;
    chebfunc_   = func;
    coeffs_.reserve(order);
    angles_.reserve(order);
    nodesTk_.reserve(order);
    cmat_ = nullptr;

    // setup some data
    setup();
}

template <class MatrixType>
void ChebyshevApproximation<MatrixType>::setup()
{
    // compute interpolation points in extents_ [a, b] (mapped from [-1,1])
    computeInterpolationPoints();
}

template <class MatrixType>
void ChebyshevApproximation<MatrixType>::buildChebyshevNodes(
    const double a, const double b, MatrixType& H)
{
    assert(order_ > 0);

    build_nodes_tm_.start();

    const int m = H.m();

    // Shift and scale matrix H so that spectrum is in range [-1, 1]
    // First make a copy of H
    MatrixType scaled_H(H);
    scaleMatrixToChebyshevInterval(a, b, scaled_H);

    // prepare to build polynomials
    const int n = order_;
    // reset nodesTk
    nodesTk_.clear();
    nodesTk_.resize(n, MatrixType("Tk", m, m));
    // Set node T_0 to identity
    nodesTk_[0].identity();
    // set node T_1 = scaled_H
    nodesTk_[1] = scaled_H;

    // define pointers for three-term recurrence
    typename std::vector<MatrixType>::iterator nodeT0, nodeT1, nodeT;
    nodeT0 = nodesTk_.begin();
    nodeT1 = nodesTk_.begin() + 1;
    nodeT  = nodesTk_.begin() + 2;

    const double m_one = -1.;
    const double two   = 2.;

    while (nodeT != nodesTk_.end())
    {
        (*nodeT) = (*nodeT0);
        (*nodeT).gemm('N', 'N', two, scaled_H, (*nodeT1), m_one);

        nodeT0 = nodeT1;
        nodeT1 = nodeT;
        nodeT++;
    }

    // done ...

    build_nodes_tm_.stop();
}

// Shift and scale matrix H so that spectrum is in range [-1, 1].
// NOTE: Matrix is modified on return
template <class MatrixType>
void ChebyshevApproximation<MatrixType>::scaleMatrixToChebyshevInterval(
    const double a, const double b, MatrixType& H)
{
    // Shift and scale matrix H so that spectrum is in range [-1, 1]
    const double scal = 2 / (b - a);
    const double shft = (b + a) / 2;
    // H = (H-shft*I)*scal
    const int m = H.m();
    std::vector<double> idmat(m, shft);
    std::vector<double> dmat(m);
    // get diagonal of H
    H.getDiagonalValues(dmat.data());
    // shift diagonal
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        m, -1., idmat.data(), dmat.data());
    // set shifted diagonal
    H.setDiagonal(dmat);
    // scale matrix
    H.scal(scal);
}

// Compute Chebyshev approximation given vectors of coefficients and Chebyshev
// nodes (polynomials)
template <class MatrixType>
MatrixType ChebyshevApproximation<MatrixType>::computeChebyshevApproximation()
{
    compute_tm_.start();

    assert(!nodesTk_.empty());
    assert(static_cast<int>(coeffs_.size()) == order_);

    // set initial approximation matrix to first polynomial approx.
    MatrixType matX(nodesTk_[0]);
    matX.scal(coeffs_[0]);

    // loop over polynomials and scale with coeffs and update matrix
    for (int k = 1; k < order_; k++)
    {
        matX.axpy(coeffs_[k], nodesTk_[k]);
    }

    compute_tm_.stop();

    return matX;
}

// Compute Chebyshev approximation given extents of approximation and matrix
// argument (for computing polynomials)
template <class MatrixType>
MatrixType ChebyshevApproximation<MatrixType>::computeChebyshevApproximation(
    const MatrixType& H, const bool recompute_coeffs)
{
    compute2_tm_.start();

    // compute Chebyshev coefficients
    if (recompute_coeffs)
    {
        computeChebyshevCoeffs();
    }
    // Shift and scale matrix H so that spectrum is in range [-1, 1]
    // First make a copy of matrix
    MatrixType scaled_H(H);
    scaleMatrixToChebyshevInterval(extents_[0], extents_[1], scaled_H);

    // define array to hold polynomials for recurrence relation
    const int m = H.m();
    std::vector<MatrixType> nodesTk(3, MatrixType("Tk", m, m));
    // Set node T_0 to identity
    nodesTk[0].identity();
    // set node T_1 = scaled_H
    nodesTk[1] = scaled_H;

    // BEGIN:
    const double m_one = -1.;
    const double two   = 2.;
    // Compute first 2 terms here:
    // 0. set initial approximation matrix to first polynomial approx.
    MatrixType mat(nodesTk[0]);
    mat.scal(coeffs_[0]);
    // 1. update with second term of approximation
    mat.axpy(coeffs_[1], nodesTk[1]);

    // define pointers for three-term recurrence
    typename std::vector<MatrixType>::iterator nodeT0, nodeT1, nodeT, begin;
    // initialize
    begin  = nodesTk.begin();
    nodeT1 = begin;
    nodeT  = begin + 1;

    // Loop to compute polynomials and update matrix
    // use pointer swapping to store polynomial updates
    for (int k = 2; k < order_; k++)
    {
        assert(coeffs_[k] == coeffs_[k]);

        // reset polynomial data pointers
        nodeT0  = nodeT1;
        nodeT1  = nodeT;
        int pos = k % 3;
        nodeT   = begin + pos;

        // Compute new polynomial node and update matrix
        (*nodeT) = (*nodeT0);
        (*nodeT).gemm('N', 'N', two, scaled_H, (*nodeT1), m_one);
        mat.axpy(coeffs_[k], (*nodeT));
    }

    // delete nodesTk;

    compute2_tm_.stop();

    return mat;
}

#ifdef MGMOL_USE_SCALAPACK
template class ChebyshevApproximation<dist_matrix::DistMatrix<DISTMATDTYPE>>;
#endif
template class ChebyshevApproximation<ReplicatedMatrix>;
