// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

#include "AndersonMix.h"
#include "Solution.h"

double residual(vector<double> diag_op, Solution& x, Solution& r)
{
    for (int i = 0; i < (int)diag_op.size(); i++)
    {
        r.setVal(i, diag_op[i] * x[i]);
    }
    double lambda = x.dotProduct(r) / x.dotProduct(x);
    cout << "lambda=" << lambda << endl;
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
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    cout << "n=" << n << endl;
    cout << "m=" << m << endl;

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

    vector<double> diag(n);
    double delta = 1. / (double)(n - 1);
    for (int i = 0; i < n; i++)
        diag[i] = i * delta;

    const double invMaxLambda = 1. / diag[n - 1];

    AndersonMix<Solution> andmix(m, beta, x);
    double normF = 1.;
    int it       = 0;

    cout << "Initial Solution:" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << x[i];
        if (i < n - 1)
            cout << ", ";
        else
            cout << endl;
    }
    do
    {
        it++;
        x.normalize();

        cout << "Trial Solution:" << endl;
        for (int i = 0; i < n; i++)
        {
            cout << x[i];
            if (i < n - 1)
                cout << ", ";
            else
                cout << endl;
        }
        normX = x.norm();
        cout << "||X||=" << normX << endl;
        x.setInvS(1. / (normX * normX));
        normF = residual(diag, x, r); // improve x
        cout << "||F||=" << normF << endl;
        r.scal(invMaxLambda);
        andmix.update(r, work, cout, true);
    } while (normF > tol && it < maxit);
    if (normF <= tol)
        cout << "Converged in " << it << " iterations" << endl;
    else
        cout << "Not Converged after " << it << " iterations" << endl;
    cout << "Solution:" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << x[i];
        if (i < n - 1)
            cout << ", ";
        else
            cout << endl;
    }

    return 0;
}
