// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef RANDOM_H
#define RANDOM_H

#include <vector>

/* Generate random numbers between (a, b] */
double generate_rand(const int a, const int b);

/* Generate random numbers between (0, 1] */
double generate_rand_num();

/* Generate a vector of random numbers between (a, b] */
std::vector<double> generate_rand(const int n);
std::vector<double> generate_rand(const int n, const int a, const int b);

template <typename DataType>
void generateRandomData(
    std::vector<DataType>& data, const DataType a, const DataType b);
#endif
