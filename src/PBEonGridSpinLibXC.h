// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_PBEONGRIDSPINLIBXC_H
#define MGMOL_PBEONGRIDSPINLIBXC_H

#ifdef USE_LIBXC

#include "Rho.h"
#include "XConGrid.h"

#include <xc.h>

#include <vector>

class Potentials;

template <class T>
class PBEonGridSpinLibXC : public XConGrid
{
    int np_;
    int myspin_;

    xc_func_type xfunc_;
    xc_func_type cfunc_;

    std::vector<double> vxc_;
    std::vector<double> exc_;
    std::vector<double> vsigma_;

    Rho<T>& rho_;

    Potentials& pot_;

public:
    PBEonGridSpinLibXC(Rho<T>& rho, Potentials& pot);

    ~PBEonGridSpinLibXC() override
    {
        xc_func_end(&xfunc_);
        xc_func_end(&cfunc_);
    }

    void update() override;

    double getExc() const override;
};
#endif

#endif
