// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ORBITALSSTEPPER_H
#define MGMOL_ORBITALSSTEPPER_H

class Ions;

template <class T>
class OrbitalsStepper
{
public:
    OrbitalsStepper() { }

    virtual ~OrbitalsStepper() { }

    virtual void setup(T&) = 0;

    virtual int updateWF(T& orbitals, Ions& ions, const double precond_factor,
        const bool orthof, T& work_orbitals, const bool accelerate,
        const bool print_res, const double atol)
        = 0;

    virtual void restartMixing() {};
};

#endif
