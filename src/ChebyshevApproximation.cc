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
#include <iostream>

Timer ChebyshevApproximation::compute_tm_("ChebyshevApproximation::compute");
Timer ChebyshevApproximation::compute2_tm_("ChebyshevApproximation::compute2");
Timer ChebyshevApproximation::build_nodes_tm_(
    "ChebyshevApproximation::build_nodes");

ChebyshevApproximation::ChebyshevApproximation(const double a, const double b,
    const int order, ChebyshevApproximationFunction* func)
{
    assert(a < b);

    extents_[0] = a;
    extents_[1] = b;
    max_order_  = order;
    order_      = order;
    chebfunc_   = func;
    coeffs_.reserve(order);
    angles_.reserve(order);
    nodesTk_ = 0;
    cmat_    = 0;

    // setup some data
    setup();
}

void ChebyshevApproximation::setup()
{
    // compute interpolation points in extents_ [a, b] (mapped from [-1,1])
    computeInterpolationPoints();
}

void ChebyshevApproximation::buildChebyshevNodes(
    const double a, const double b, dist_matrix::DistMatrix<DISTMATDTYPE>& H)
{
    assert(order_ > 0);

    build_nodes_tm_.start();

    const int m = H.m();

    // Shift and scale matrix H so that spectrum is in range [-1, 1]
    // First make a copy of H
    dist_matrix::DistMatrix<DISTMATDTYPE> scaled_H(H);
    scaleMatrixToChebyshevInterval(a, b, scaled_H);

    // reset nodesTk
    if (nodesTk_)
    {
        delete nodesTk_;
        nodesTk_ = 0;
    }

    // prepare to build polynomials
    const int n = order_;

    std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>>* nodes
        = new std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>>(
            n, dist_matrix::DistMatrix<DISTMATDTYPE>("Tk", m, m));
    // Set node T_0 to identity
    nodes->at(0).identity();
    // set node T_1 = scaled_H
    nodes->at(1) = scaled_H;

    // define pointers for three-term recurrence
    std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>>::iterator nodeT0, nodeT1,
        nodeT;
    nodeT0 = nodes->begin();
    nodeT1 = nodes->begin() + 1;
    nodeT  = nodes->begin() + 2;

    const double m_one = -1.;
    const double two   = 2.;

    while (nodeT != nodes->end())
    {
        (*nodeT) = (*nodeT0);
        (*nodeT).gemm('N', 'N', two, scaled_H, (*nodeT1), m_one);

        nodeT0 = nodeT1;
        nodeT1 = nodeT;
        nodeT++;
    }

    // done ... store pointer to nodes
    nodesTk_ = nodes;

    build_nodes_tm_.stop();
}

// Shift and scale matrix H so that spectrum is in range [-1, 1].
// NOTE: Matrix is modified on return
void ChebyshevApproximation::scaleMatrixToChebyshevInterval(
    const double a, const double b, dist_matrix::DistMatrix<DISTMATDTYPE>& H)
{
    // Shift and scale matrix H so that spectrum is in range [-1, 1]
    const double scal = 2 / (b - a);
    const double shft = (b + a) / 2;
    // H = (H-shft*I)*scal
    const int m = H.m();
    std::vector<double> idmat(m, shft);
    std::vector<double> dmat(m);
    // get diagonal of H
    H.getDiagonalValues(&dmat[0]);
    // shift diagonal
    MPaxpy(m, -1., &idmat[0], &dmat[0]);
    // set shifted diagonal
    H.setDiagonal(dmat);
    // scale matrix
    H.scal(scal);
}

// Compute Chebyshev approximation given vectors of coefficients and Chebyshev
// nodes (polynomials)
dist_matrix::DistMatrix<DISTMATDTYPE>
ChebyshevApproximation::computeChebyshevApproximation()
{
    compute_tm_.start();

    assert(nodesTk_ != 0);
    assert((int)coeffs_.size() == order_);
    //    assert((int)nodesTk_->size() == order_);

    // set initial approximation matrix to first polynomial approx.
    dist_matrix::DistMatrix<DISTMATDTYPE> matX(nodesTk_->at(0));
    matX.scal(coeffs_[0]);

    // loop over polynomials and scale with coeffs and update matrix
    for (int k = 1; k < order_; k++)
    {
        matX.axpy(coeffs_[k], nodesTk_->at(k));
    }

    compute_tm_.stop();

    return matX;
}

// Compute Chebyshev approximation given extents of approximation and matrix
// argument (for computing polynomials)
dist_matrix::DistMatrix<DISTMATDTYPE>
ChebyshevApproximation::computeChebyshevApproximation(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& H, const bool recompute_coeffs)
{
    compute2_tm_.start();

    // compute Chebyshev coefficients
    if (recompute_coeffs)
    {
        computeChebyshevCoeffs();
    }
    // Shift and scale matrix H so that spectrum is in range [-1, 1]
    // First make a copy of matrix
    dist_matrix::DistMatrix<DISTMATDTYPE> scaled_H(H);
    scaleMatrixToChebyshevInterval(extents_[0], extents_[1], scaled_H);

    // define array to hold polynomials for recurrence relation
    const int m = H.m();
    std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>>* nodesTk
        = new std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>>(
            3, dist_matrix::DistMatrix<DISTMATDTYPE>("Tk", m, m));
    // Set node T_0 to identity
    nodesTk->at(0).identity();
    // set node T_1 = scaled_H
    nodesTk->at(1) = scaled_H;

    // BEGIN:
    const double m_one = -1.;
    const double two   = 2.;
    // Compute first 2 terms here:
    // 0. set initial approximation matrix to first polynomial approx.
    dist_matrix::DistMatrix<DISTMATDTYPE> mat(nodesTk->at(0));
    mat.scal(coeffs_[0]);
    // 1. update with second term of approximation
    mat.axpy(coeffs_[1], nodesTk->at(1));

    // define pointers for three-term recurrence
    std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>>::iterator nodeT0, nodeT1,
        nodeT, begin;
    // initialize
    begin  = nodesTk->begin();
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

    delete nodesTk;

    compute2_tm_.stop();

    return mat;
}
