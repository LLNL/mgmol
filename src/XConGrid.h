// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_XCONGRID_H
#define MGMOL_XCONGRID_H

#include "Timer.h"

class XConGrid
{

public:
    static Timer get_xc_tm_;

    XConGrid(){};

    virtual ~XConGrid(){};

    virtual void update() = 0;

    virtual double getExc() const = 0;
};

#endif
