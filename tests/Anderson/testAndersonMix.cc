// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cmath>
#include <iostream>
#include <vector>

#include "AndersonMix.h"
#include "Solution.h"

double residual(std::vector<double> diag_op, Solution& x, Solution& r)
{
    for (int i = 0; i < (int)diag_op.size(); i++)
    {
        r.setVal(i, diag_op[i] * x[i]);
    }
    double lambda = x.dotProduct(r) / x.dotProduct(x);
    std::cout << "lambda=" << lambda << std::endl;
    for (int i = 0; i < (int)diag_op.size(); i++)
    {
        r.setVal(i, lambda * x[i] - r[i]);
    }

    double normR = r.dotProduct(r);

    return normR;
}

// TODO this test does not check anything
int main(int argc, char** argv)
{
    if (argc < 3)
    {
        std::cerr << "Insufficient number of arguments" << std::endl;
    }
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    std::cout << "n=" << n << std::endl;
    std::cout << "m=" << m << std::endl;

    const int maxit  = 100;
    const double tol = 1.e-8;

    double beta = 1.;

    Solution x(n);
    x.initRand(n);
    x.setVal(0, 1.);
    x.setVal(1, 2.);

    Solution work(n);

    double normX = sqrt(x.dotProduct(x));
    x.scal(1. / normX);

    Solution r(n);

    std::vector<double> diag(n);
    double delta = 1. / (double)(n - 1);
    for (int i = 0; i < n; i++)
        diag[i] = i * delta;

    const double invMaxLambda = 1. / diag[n - 1];

    AndersonMix<Solution> andmix(m, beta, x);
    double normF = 1.;
    int it       = 0;

    std::cout << "Initial Solution:" << std::endl;
    for (int i = 0; i < n; i++)
    {
        std::cout << x[i];
        if (i < n - 1)
            std::cout << ", ";
        else
            std::cout << std::endl;
    }
    do
    {
        it++;
        x.normalize();

        std::cout << "Trial Solution:" << std::endl;
        for (int i = 0; i < n; i++)
        {
            std::cout << x[i];
            if (i < n - 1)
                std::cout << ", ";
            else
                std::cout << std::endl;
        }
        normX = x.norm();
        std::cout << "||X||=" << normX << std::endl;
        x.setInvS(1. / (normX * normX));
        normF = residual(diag, x, r); // improve x
        std::cout << "||F||=" << normF << std::endl;
        r.scal(invMaxLambda);
        andmix.update(r, work, std::cout, true);
    } while (normF > tol && it < maxit);
    if (normF <= tol)
        std::cout << "Converged in " << it << " iterations" << std::endl;
    else
        std::cout << "Not Converged after " << it << " iterations" << std::endl;
    std::cout << "Solution:" << std::endl;
    for (int i = 0; i < n; i++)
    {
        std::cout << x[i];
        if (i < n - 1)
            std::cout << ", ";
        else
            std::cout << std::endl;
    }

    return 0;
}
