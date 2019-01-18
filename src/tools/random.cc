// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "random.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

/* Generate random numbers between (a, b] */
double generate_rand(const int a, const int b)
{
    double val = (b - a) * (double)rand() / (double)RAND_MAX + a;
    return val;
}

/* Generate random numbers between (0, 1] */
double generate_rand_num()
{
    double val = (double)rand() / (double)RAND_MAX;
    return val;
}

/* Generate a vector of random numbers between (a, b] */
std::vector<double> generate_rand(const int n, const int a, const int b)
{
    std::vector<double> vec(n);
    for (int i = 0; i < n; i++)
        vec[i] = generate_rand(a, b);

    return vec;
}
std::vector<double> generate_rand(const int n)
{
    std::vector<double> vec(n);
    std::generate(vec.begin(), vec.end(), &generate_rand_num);

    return vec;
}
