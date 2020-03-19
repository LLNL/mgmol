// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <iostream>
#include <stdlib.h>
#include <vector>

void splint(
    double* xa, double* ya, double* y2a, const int n, double x, double* y)
{
    int klo = 0;
    int khi = n - 1;

    while (khi - klo > 1)
    {
        int k = (khi + klo) >> 1;
        if (xa[k] > x)
            khi = k;
        else
            klo = k;
    }
    double h = xa[khi] - xa[klo];
    if (h == 0.)
    {
        std::cout << "Bad xa input" << std::endl;
        exit(0);
    }
    double a = (xa[khi] - x) / h;
    double b = (x - xa[klo]) / h;
    *y       = a * ya[klo] + b * ya[khi]
         + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h)
               / 6.;
}

