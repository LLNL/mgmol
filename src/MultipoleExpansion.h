// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef MULTIPOLEEXPANSION_H
#define MULTIPOLEEXPANSION_H

#include "GridFunc.h"
#include "MPIdata.h"
#include "Rho.h"
#include "Vector3D.h"

#include <vector>

// template <typename T>
class MultipoleExpansion
{
private:
    const pb::Grid& grid_;
    Vector3D origin_;
    Vector3D origin_cell_;
    Vector3D cell_;
    short bc_[3];

    short order_;
    double qtotal_;
    Vector3D dipole_moment_;
    double quadrupole_moment_[6];

    std::vector<short> space_dim_;

    bool onpe0_;

    void get_dipole(RHODTYPE* rho);
    void get_monopole(RHODTYPE* rho);

    void get_dipole(const pb::GridFunc<RHODTYPE>& rho);
    void get_monopole(const pb::GridFunc<RHODTYPE>& rho);
    void computeQuadrupole(const pb::GridFunc<RHODTYPE>& rho);
    void resetOriginToChargeCenter(const pb::GridFunc<RHODTYPE>& rho);

public:
    // constructor
    MultipoleExpansion(const pb::Grid& mygrid, const short bc[3],
        const Vector3D&, const Vector3D&);

    void setup(RHODTYPE* rho);
    void setup(pb::GridFunc<RHODTYPE>& rho);
    template <typename T>
    void expand(pb::GridFunc<T>& func);
    template <typename T>
    void expand2d(pb::GridFunc<T>& func);
    template <typename T>
    void expand3d(pb::GridFunc<T>& func);

    void setOrder(const short order) { order_ = order; }
};

#endif
