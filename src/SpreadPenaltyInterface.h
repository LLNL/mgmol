// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef SpreadPenaltyInterface_H
#define SpreadPenaltyInterface_H

class LocGridOrbitals;

class SpreadPenaltyInterface
{
public:
    SpreadPenaltyInterface() {}

    virtual ~SpreadPenaltyInterface() {}

    // add penalty functional contribution to residual
    virtual void addResidual(LocGridOrbitals& phi, LocGridOrbitals& res) = 0;
    virtual double evaluateEnergy(const LocGridOrbitals& phi)            = 0;
};

#endif
