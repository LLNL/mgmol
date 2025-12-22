// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_SpreadPenaltyInterface_H
#define MGMOL_SpreadPenaltyInterface_H

template <class T>
class SpreadPenaltyInterface
{
public:
    SpreadPenaltyInterface() { }

    virtual ~SpreadPenaltyInterface() { }

    // add penalty functional contribution to residual
    virtual void addResidual(T& phi, T& res)    = 0;
    virtual double evaluateEnergy(const T& phi) = 0;
};

#endif
