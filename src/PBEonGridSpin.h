// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_PBEONGRIDSPIN_H
#define MGMOL_PBEONGRIDSPIN_H

#include "Mesh.h"
#include "Rho.h"
#include "XConGrid.h"

#include "PBEFunctional.h"

#include <vector>

class Potentials;

template <class T>
class PBEonGridSpin : public XConGrid
{
    int np_;
    int myspin_;
    std::vector<double> vxc_;
    PBEFunctional* pbe_;
    Rho<T>& rho_;

    Potentials& pot_;

public:
    PBEonGridSpin(Rho<T>& rho, Potentials& pot);

    ~PBEonGridSpin() override { delete pbe_; }

    void update() override;

    double getExc() const override;
};

#endif
