// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <iostream>
#include <vector>

void spline(
    double* x, double* y, const int n, double yp1, double ypn, double* y2)
{
    double p, qn, sig, un;

    vector<double> u(n);
    if (yp1 > 0.99e30)
    {
        cout << "Set ypp to 0 on left" << endl;
        y2[0] = u[0] = 0.;
    }
    else
    {
        y2[0] = -0.5;
        u[0]  = (3. / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }
    for (int i = 1; i < n - 1; i++)
    {
        sig   = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p     = sig * y2[i - 1] + 2.;
        y2[i] = (sig - 1.) / p;
        u[i]  = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
               - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (6. * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    if (ypn > 0.99e30)
        qn = un = 0.0;
    else
    {
        qn = 0.5;
        un = (3. / (x[n - 1] - x[n - 2]))
             * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
    }
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.);
    for (int k = n - 2; k >= 0; k--)
        y2[k] = y2[k] * y2[k + 1] + u[k];

    std::cout << "y2[0]=" << y2[0] << std::endl;
}
