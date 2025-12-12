// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_RADIALINTER_H
#define MGMOL_RADIALINTER_H

#include "RadialMeshFunction.h"

class RadialInter : public RadialMeshFunction
{
private:
public:
    RadialInter(const std::vector<double>& x) : RadialMeshFunction(x) { }

    RadialInter() : RadialMeshFunction() { }

    RadialInter(std::vector<double>& x, std::vector<std::vector<double>>& y)
        : RadialMeshFunction(x, y)
    {
    }

    ~RadialInter() override { }

    double linint(const double x, const int j = 0) const;
    double cubint(const double x, const int j = 0) const;
};

#endif
