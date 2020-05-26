// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_CHEBYSHEV_APPROX_FUNC_H
#define MGMOL_CHEBYSHEV_APPROX_FUNC_H

#include <vector>

class ChebyshevApproximationFunction
{

private:
public:
    ChebyshevApproximationFunction(){};
    ~ChebyshevApproximationFunction(){}; // destructor

    virtual std::vector<double> eval(const std::vector<double>& x) = 0;
};

#endif
