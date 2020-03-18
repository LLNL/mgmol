// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "BlacsContext.h"
#include "DistMatrix.h"
#include "MGmol_MPI.h"

#include "catch.hpp"
#include <iostream>

TEST_CASE("Check condition number of DistMatrix using Cholesky decomposition",
    "[condition_dist_matrix]")
{
    MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

    int nprow = 2;
    int npcol = 2;
    dist_matrix::BlacsContext bc(MPI_COMM_WORLD, nprow, npcol);

    const int n        = 51;
    const double shift = 0.2;

    // block size
    const int nb = 7;

    // construct diagonal 10x10 matrix
    dist_matrix::DistMatrix<double> S("S", bc, n, n, nb, nb);
    for (int i = 0; i < n; i++)
        S.setVal(i, i, 0.1 * i + shift);

    const double expected_cond = (0.1 * (n - 1) + shift) / (shift);

    dist_matrix::DistMatrix<double> LS("LS", bc, n, n, nb, nb);

    // Cholesky decomposition of S
    LS       = S;
    int info = LS.potrf('l');
    REQUIRE(info == 0);

    double anorm   = S.norm('1');
    double invcond = LS.pocon('l', anorm);
    double cond    = 1. / invcond;

    const double tol = 1.e-3;
    CHECK(cond == Approx(expected_cond).epsilon(tol));
}
