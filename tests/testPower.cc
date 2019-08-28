// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "LocalVector.h"
#include "SquareLocalMatrices.h"

#include "Power.h"

#include <boost/test/unit_test.hpp>
#include <iostream>

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(power, *utf::tolerance(1e-3))
{
    // construct diagonal 10x10 matrix
    const int n        = 10;
    const double shift = -3.;

    SquareLocalMatrices<double> A(1, n);
    for (int i = 0; i < n; i++)
        A.setVal(i, i, 10. * i + shift);

    double emin;
    double emax;
    Power<LocalVector<double>, SquareLocalMatrices<double>> power(n);

    power.computeEigenInterval(A, emin, emax, 1.e-3, true);

    BOOST_TEST(emin == shift);

    const double expected_emax = 10. * (n - 1) + shift;
    BOOST_TEST(emax == expected_emax);
}
