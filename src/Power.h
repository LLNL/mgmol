// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef MGMOL_POWER_H
#define MGMOL_POWER_H

#include "Timer.h"
#include "random.h"
#include <iostream>

template <class VECTOR, class MATRIX>
class Power
{
private:
    static Timer compute_tm_;

    // use shift to target highest or lowest eigenvalue
    double shift_;

    VECTOR vec1_;
    VECTOR vec2_;

public:
    Power(const int n)
        : shift_(0.), vec1_(generate_rand(n)), vec2_(generate_rand(n))
    {
    }

    void computeEigenInterval(const MATRIX& A, double& emin, double& emax,
        const double epsilon = 1.e-2, const bool verbose = false);

    static void printTimers(std::ostream& os) { compute_tm_.print(os); }
};

template <class VECTOR, class MATRIX>
Timer Power<VECTOR, MATRIX>::compute_tm_("Power::compute");

#endif
