// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ORBITALSEXTRAPOLATIONORDER2_H
#define MGMOL_ORBITALSEXTRAPOLATIONORDER2_H

#include "OrbitalsExtrapolation.h"

class OrbitalsExtrapolationOrder2 : public OrbitalsExtrapolation
{
public:
    
    OrbitalsExtrapolationOrder2()
    {
        extrapolated_H_=false;
    }
    
    void extrapolate_orbitals(LocGridOrbitals** orbitals, 
                              LocGridOrbitals* new_orbitals);

    bool extrapolatedH()const
    {
        return extrapolated_H_;
    }
private:
    
    bool extrapolated_H_;

};

#endif
