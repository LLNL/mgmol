// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_PBEONGRIDLIBXC_H
#define MGMOL_PBEONGRIDLIBXC_H

#ifdef USE_LIBXC

#include "Rho.h"
#include "XConGrid.h"

#include <xc.h>

#include <vector>

class Potentials;

template <class T>
class PBEonGridLibXC : public XConGrid
{
    int np_;
    xc_func_type xfunc_;
    xc_func_type cfunc_;
    std::vector<double> exc_;
    std::vector<double> vxc_;
    std::vector<double> vsigma_;
    Rho<T>& rho_;

    Potentials& pot_;

public:
    PBEonGridLibXC(Rho<T>& rho, Potentials& pot) : rho_(rho), pot_(pot)
    {
        np_         = rho.rho_[0].size();
        int func_id = XC_GGA_X_PBE;
        if (xc_func_init(&xfunc_, func_id, XC_UNPOLARIZED) != 0)
        {
            std::cerr << "Functional " << func_id << " not found" << std::endl;
        }
        func_id = XC_GGA_C_PBE;
        if (xc_func_init(&cfunc_, func_id, XC_UNPOLARIZED) != 0)
        {
            std::cerr << "Functional " << func_id << " not found" << std::endl;
        }
        exc_.resize(np_);
        vxc_.resize(np_);
        vsigma_.resize(np_);
    }

    ~PBEonGridLibXC() override
    {
        xc_func_end(&xfunc_);
        xc_func_end(&cfunc_);
    }

    void update() override;

    double getExc() const override;
};

#endif

#endif
