// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Poisson.h 1304 2015-09-03 00:31:15Z osei $
#ifndef POISSONINTERFACE_H_
#define POISSONINTERFACE_H_

#include "Timer.h"

#include <iostream>

class PoissonInterface
{
protected:
    static Timer poisson_tm_;

public:
    virtual ~PoissonInterface() { }
    static void printTimers(std::ostream& os) { poisson_tm_.print(os); }
};

#endif
