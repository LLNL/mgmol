//
// test Replicated2DistMatrix
//
#include "BlacsContext.h"
#include "DistMatrix.h"
#include "MGmol_MPI.h"
#include "SquareLocalMatrices.h"

#include "catch.hpp"

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

TEST_CASE(
    "Check that distributed matrices can be build from replicated matrices",
    "[replicated_to_dist_matrix]")
{
    int mype;
    int npes;

    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    INFO("This example to set up to use only 4 processes.");
    REQUIRE(npes == 4);
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    int nprow = 2;
    int npcol = 2;

    const int n  = 11;
    const int nb = 3;

    int m    = n;
    int mb_a = nb;
    int nb_a = nb;

    MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

    dist_matrix::BlacsContext bc(MPI_COMM_WORLD, nprow, npcol);

    dist_matrix::DistMatrix<double> distm("a", bc, m, n, mb_a, nb_a);

    if (mype == 0)
    {
        std::cout << " m x n / mb x nb = " << distm.m() << "x" << distm.n()
                  << " / " << distm.mb() << "x" << distm.nb() << " / "
                  << std::endl;
    }

    // setup an nxn local matrix
    SquareLocalMatrices<double, MemorySpace::Host> replicated(1, n);

    std::vector<double> tmp(n * m);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            tmp[i + j * m] = 10. * (i + 1) + j + 1;

    replicated.setValues(tmp.data(), m);

    // distribute replicated matrix
    distm.initFromReplicated(replicated.getRawPtr(), n);

    // convert back to a replicated matrix
    SquareLocalMatrices<double, MemorySpace::Host> result(1, n);
    distm.allgather(result.getRawPtr(), n);

    if (mype == 0)
    {
        std::cout << std::setprecision(2);
        result.print(std::cout);
    }

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
        {
            double valbefore = replicated.getVal(i, j);
            double valafter  = result.getVal(i, j);
            CHECK(valbefore == Approx(valafter).epsilon(1e-8));
        }
}
