// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ORTHOANDERSON_MIX
#define MGMOL_ORTHOANDERSON_MIX

#include "AndersonMix.h"

template <class T>
class OrthoAndersonMix : public AndersonMix<T>
{
private:
    const int m_;
    T& x_; // current trial solution

protected:
    void postprocessUpdate() override;

public:
    OrthoAndersonMix(const int m, const double beta, T& x)
        : AndersonMix<T>(m, beta, x), m_(m), x_(x)
    {
    }

    ~OrthoAndersonMix() override {};
};

#endif
