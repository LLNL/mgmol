// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef ORBITALSSTEPPER_H
#define ORBITALSSTEPPER_H

class OrbitalsStepper
{
public:

    OrbitalsStepper()
    {}
    
    virtual ~OrbitalsStepper()
    {}
    
    virtual void setup(LocGridOrbitals&)=0;
    
    virtual int update(LocGridOrbitals& orbitals,
        Ions& ions,
        const double precond_factor,
        const bool orthof,
        LocGridOrbitals& work_orbitals,
        const bool accelerate,
        const bool print_res,
        const double atol)=0;
};

#endif  
