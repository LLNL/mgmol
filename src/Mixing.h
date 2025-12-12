// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_MIXING
#define MGMOL_MIXING

#include <iostream>

template <class T>
class Mixing
{
public:
    Mixing() {};

    virtual ~Mixing() {};

    virtual void update(T& res, T& work, std::ostream& os, const bool verbose)
        = 0;
    virtual void restart(void) = 0;
};

#endif
