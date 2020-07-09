// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "random.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

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

template <typename DataType>
void generateRandomData(
    std::vector<DataType>& data, const DataType minv, const DataType maxv)
{
    typedef boost::minstd_rand rng_type;
    typedef boost::uniform_real<> distribution_type;

    int seed = 113;
    rng_type rng(seed);
    distribution_type nd(minv, maxv);
    boost::variate_generator<rng_type, distribution_type> gen(rng, nd);

    for (auto& d : data)
        d = gen();
}

template void generateRandomData(
    std::vector<double>& data, const double minv, const double maxv);
template void generateRandomData(
    std::vector<float>& data, const float minv, const float maxv);
