// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef MGMOL_XYZOPS_H
#define MGMOL_XYZOPS_H

#include "Timer.h"

#include <vector>

template <class OrbitalsType>
class XYZOps
{
private:
    static Timer compute_tm_;

public:
    static void compute(
        const OrbitalsType& orbitals, std::vector<std::vector<double>>& a);

    static void printTimers(std::ostream& os) { compute_tm_.print(os); }
};

template <class OrbitalsType>
Timer XYZOps<OrbitalsType>::compute_tm_("SinCosOps::compute_tm");

#endif
