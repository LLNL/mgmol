//
// test DistMatrix
//
// multiply a matrix a(m,k) by b(k,n) to get c(m,n)
// using blocks of size (mb,nb) on a process grid (nprow,npcol)
//
// use: testDistMatrix [-check] input_file
// input:
// nprow npcol
// m_a n_a mb_a nb_a transa
// m_b n_b mb_b nb_b transb
// m_c n_c mb_c nb_c
//
#include "BlacsContext.h"
#include "DistMatrix.h"
#include "DistVector.h"
#include "MGmol_MPI.h"

#include "catch.hpp"

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

TEST_CASE("Check DistVector", "[dist_vector]")
{
    int mype;
    int npes;

    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    INFO("This example to set up to use only 4 processes");
    REQUIRE(npes == 4);
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    int nprow = 2;
    int npcol = 2;

    const int m  = 21;
    const int nb = 6;

    {
        MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

        dist_matrix::BlacsContext bc(MPI_COMM_WORLD, nprow, npcol);
        dist_matrix::DistMatrix<DISTMATDTYPE>::setBlockSize(nb);

        dist_matrix::DistMatrix<DISTMATDTYPE>::setDefaultBlacsContext(&bc);

        dist_matrix::DistVector<double> distv("v", m);

        std::vector<double> stdv;
        stdv.reserve(m);
        for (int i = 0; i < m; i++)
            stdv.push_back(1.);

        distv.assign(stdv);

        double norm2          = distv.nrm2();
        double expected_norm2 = std::sqrt((double)m);

        CHECK(norm2 == Approx(expected_norm2).epsilon(0.000001));

        distv.normalize();

        double normn = distv.nrm2();
        CHECK(normn == Approx(1.).epsilon(0.000001));

        dist_matrix::DistVector<double> distv2(stdv);

        double diff2 = distv2.scaledDiff2(distv, norm2);
        CHECK(diff2 == Approx(0.).margin(0.000001));
    }
}
