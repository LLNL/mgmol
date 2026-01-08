// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "DensityMatrix.h"
#include "GramMatrix.h"

#ifdef MGMOL_USE_SCALAPACK
#include "BlacsContext.h"
#include "DistMatrix.h"
#else
#include "ReplicatedMatrix.h"
#endif

#include "catch.hpp"

#include <vector>

TEST_CASE(
    "Check functionalities of class DensityMatrix", "[functions_DensityMatrix")
{
#ifdef MGMOL_USE_SCALAPACK
    typedef dist_matrix::DistMatrix<double> MatrixType;
#else
    typedef ReplicatedMatrix MatrixType;
#endif
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int npes;
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

#ifdef MGMOL_USE_SCALAPACK
    INFO("This example to set up to use only 4 processes");
    REQUIRE(npes == 4);

    int nprow = 2;
    int npcol = 2;
    dist_matrix::BlacsContext bc(MPI_COMM_WORLD, nprow, npcol);
    dist_matrix::DistMatrix<DISTMATDTYPE>::setDefaultBlacsContext(&bc);
#endif
    const int n = 223;

    std::vector<double> raw_mat(n * n, 0.);

    raw_mat[0] = 2.;
    raw_mat[1] = -1.;
    for (int i = 1; i < n - 1; i++)
    {
        raw_mat[(n + 1) * i]     = 2.;
        raw_mat[(n + 1) * i + 1] = -1.;
        raw_mat[(n + 1) * i - 1] = -1.;
    }
    raw_mat[(n + 1) * (n - 1)]     = 2.;
    raw_mat[(n + 1) * (n - 1) - 1] = -1.;

    MatrixType matK("K", n);
    matK.init(raw_mat.data(), n);
    if (myrank == 0) std::cout << "K" << std::endl;
    matK.print(std::cout, 0, 0, 5, 5);

    const double two_third = 2. / 3.;
    const double one_sixth = 1. / 6.;
    raw_mat[0]             = two_third;
    raw_mat[1]             = one_sixth;
    for (int i = 1; i < n - 1; i++)
    {
        raw_mat[(n + 1) * i]     = two_third;
        raw_mat[(n + 1) * i + 1] = one_sixth;
        raw_mat[(n + 1) * i - 1] = one_sixth;
    }
    raw_mat[(n + 1) * (n - 1)]     = 2.;
    raw_mat[(n + 1) * (n - 1) - 1] = one_sixth;

    MatrixType matM("M", n);
    matM.init(raw_mat.data(), n);
    if (myrank == 0) std::cout << "M" << std::endl;
    matM.print(std::cout, 0, 0, 5, 5);

    // new Gram matrix
    GramMatrix<MatrixType> gram(n);

    gram.setMatrix(matM, 0);

    // compute Cholesky decomposition
    gram.updateLS();
    const MatrixType& ls = gram.getCholeskyL();

    // setup density matrix
    DensityMatrix<MatrixType> dm(n);

    dm.setMatrix(matK);

    dm.stripS(ls);
    dm.dressUpS(ls);

    const MatrixType& newM = dm.getMatrix();
    if (myrank == 0) std::cout << "new M" << std::endl;
    newM.print(std::cout, 0, 0, 5, 5);

    matK -= newM;

    double normK = matK.norm('m');
    CHECK(normK == Approx(0.).margin(1.e-14));
}
