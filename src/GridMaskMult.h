// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef GRIDMASKMULT_H
#define GRIDMASKMULT_H

#include "GridMask.h"

class GridMaskMult : public GridMask
{
    template <typename T>
    void applyPrivate(
        T* u, const unsigned short level, const unsigned short iloc);
    template <typename T>
    void applyPrivate(pb::GridFunc<T>& gu, const unsigned short level,
        const unsigned short iloc);

public:
    GridMaskMult(const unsigned short nclevels, const unsigned short subdivx,
        const pb::Grid& mygrid)
        : GridMask(nclevels, subdivx, mygrid) {};

    void apply(float* u, const unsigned short level, const unsigned short iloc,
        const bool /*first_application*/ = false) override
    {
        applyPrivate(u, level, iloc);
    }

    void apply(double* u, const unsigned short level, const unsigned short iloc,
        const bool /*first_application*/ = false) override
    {
        applyPrivate(u, level, iloc);
    }

    void apply(pb::GridFunc<double>& gu, const unsigned short level,
        const unsigned short iloc) override
    {
        applyPrivate(gu, level, iloc);
    }

    void apply(pb::GridFunc<float>& gu, const unsigned short level,
        const unsigned short iloc) override
    {
        applyPrivate(gu, level, iloc);
    }
};

#endif
