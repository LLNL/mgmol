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

#include <iostream>

int main(int argc, char** argv)
{
    std::cout << "Test Power" << std::endl;

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

    return 0;
}
