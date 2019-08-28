// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "BlacsContext.h"
#include "DistMatrix.h"
#include "DistVector.h"
#include "MGmol_MPI.h"
#include "Power.h"

#include <boost/test/unit_test.hpp>
#include <iostream>

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(power_dist_matrix, *utf::tolerance(1e-3))
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

    int nprow = 2;
    int npcol = 2;
    dist_matrix::BlacsContext bc(MPI_COMM_WORLD, nprow, npcol);

    const int n        = 10;
    const double shift = -3.;

    // block size
    const int nb = 6;

    // construct diagonal 10x10 matrix
    dist_matrix::DistMatrix<double> A("A", bc, n, n, nb, nb);
    for (int i = 0; i < n; i++)
        A.setVal(i, i, 10. * i + shift);

    // DistVector will be generated in Power
    // and need default values
    dist_matrix::DistMatrix<DISTMATDTYPE>::setDefaultBlacsContext(&bc);
    double emin;
    double emax;
    Power<dist_matrix::DistVector<double>, dist_matrix::DistMatrix<double>>
        power(n);

    power.computeEigenInterval(A, emin, emax, 1.e-3, (myrank == 0));

    BOOST_TEST(emin == shift);

    const double expected_emax = 10. * (n - 1) + shift;
    BOOST_TEST(emax == expected_emax);
}
