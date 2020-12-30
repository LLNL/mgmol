// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "LocalVector.h"
#include "SquareLocalMatrices.h"

#include "Power.h"

#include "catch.hpp"
#include <iostream>

TEST_CASE("Compute eigenvalues using the power method", "[power]")
{
    // construct diagonal 10x10 matrix
    const int n        = 10;
    const double shift = -3.;

    SquareLocalMatrices<double, MemorySpace::Host> A(1, n);
    std::vector<double> tmp(n * n);
    memset(tmp.data(), 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++)
        tmp[i + n * i] = 10. * i + shift;

    A.setValues(tmp.data(), n);

    double emin;
    double emax;
    Power<LocalVector<double, MemorySpace::Host>,
        SquareLocalMatrices<double, MemorySpace::Host>>
        power(n);

    power.computeEigenInterval(A, emin, emax, 1.e-3, true);

    CHECK(emin == Approx(shift).epsilon(1e-3));

    const double expected_emax = 10. * (n - 1) + shift;
    CHECK(emax == Approx(expected_emax).epsilon(1e-3));
}
