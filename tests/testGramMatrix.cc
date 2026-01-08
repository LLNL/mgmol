// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "GramMatrix.h"
#include "ReplicatedMatrix.h"

#ifdef MGMOL_USE_SCALAPACK
#include "BlacsContext.h"
#include "DistMatrix.h"
#endif

#include "catch.hpp"

#include <vector>

TEST_CASE("Check functionalities of class GramMatrix", "[functions_GramMatrix")
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
    const int n = 213;

    // build matrix with diagonal elements equal to 1.
    // and off diagonal elements such that its eigenvalues
    // are guaranteed to be between 0 and 2
    MatrixType matA("A", n);

    matA.setRandom(-1. / (double)n, 1. / (double)n);

    MatrixType matB("B", n);
    matB = matA;
    matB.transpose(0.5, matA, 0.5);

    std::vector<double> diag(n, 2.);
    matB.setDiagonal(diag);

    if (myrank == 0) std::cout << "Matrix B:" << std::endl;
    matB.print(std::cout, 0, 0, 5, 5);

    double normB = matB.norm('m');
    if (myrank == 0) std::cout << "Norm B = " << normB << std::endl;

    // new Gram matrix
    GramMatrix<MatrixType> gram(n);

    // initialize Gram matrix with "matB"
    gram.setMatrix(matB, 0);

    // compute Cholesky decomposition
    gram.updateLS();

    MatrixType lmat(gram.getCholeskyL());
    lmat.trset('l');

    // compute product L*L^T
    matA.gemm('N', 'T', 1., lmat, lmat, 0.);
    if (myrank == 0) std::cout << "L*L^T:" << std::endl;
    matA.print(std::cout, 0, 0, 5, 5);

    matA -= matB;
    if (myrank == 0) std::cout << "L*L^T-B:" << std::endl;
    matA.print(std::cout, 0, 0, 5, 5);

    double normA = matA.norm('m');
    if (myrank == 0) std::cout << "Norm A-B = " << normA << std::endl;

    CHECK(normA == Approx(0.).margin(1.e-14));

    // Loewdin
    MatrixType loewdin("Loewdin", n);
    std::shared_ptr<MatrixType> invloewdin;
    invloewdin.reset(new MatrixType("InvLoewdin", n));
    gram.computeLoewdinTransform(loewdin, invloewdin, 0);
    if (myrank == 0) std::cout << "Loewdin" << std::endl;
    loewdin.print(std::cout, 0, 0, 5, 5);

    matA.gemm('N', 'N', 1., loewdin, *invloewdin, 0.);
    if (myrank == 0) std::cout << "B^-1/2*B^1/2:" << std::endl;
    matA.print(std::cout, 0, 0, 5, 5);

    MatrixType id("Id", n);
    id.identity();
    matA -= id;

    normA = matA.norm('m');
    CHECK(normA == Approx(0.).margin(1.e-14));

    matA.gemm('N', 'N', 1., *invloewdin, *invloewdin, 0.);
    if (myrank == 0) std::cout << "B^1/2*B^1/2:" << std::endl;
    matA.print(std::cout, 0, 0, 5, 5);

    matA -= matB;

    normA = matA.norm('m');
    CHECK(normA == Approx(0.).margin(1.e-12));
}
