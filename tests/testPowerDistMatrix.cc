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

#include <iostream>

int main(int argc, char** argv)
{
    int mpirc = MPI_Init(&argc, &argv);
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0) std::cout << "Test PowerDistMatrix" << std::endl;

    {
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

        const double tol = 1.e-3;
        if (std::abs(emin - shift) > tol)
        {
            std::cerr << "Lowest eigenvalue (" << emin << ")is incorrect!!!"
                      << std::endl;
            std::cerr << "Error: " << emin - shift << std::endl;
            return 1;
        }

        const double expected_emax = 10. * (n - 1) + shift;
        if (std::abs(emax - expected_emax) > tol)
        {
            std::cerr << "Highest eigenvalue (" << emax << ")is incorrect!!!"
                      << std::endl;
            std::cerr << "Error: " << emax - expected_emax << std::endl;
            return 1;
        }
    }

    mpirc = MPI_Finalize();

    // return 0 for SUCCESS
    return 0;
}
