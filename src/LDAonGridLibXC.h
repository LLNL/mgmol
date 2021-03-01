// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_LDAONGRIDLIBXC_H
#define MGMOL_LDAONGRIDLIBXC_H

#ifdef USE_LIBXC

#include "MGmol_MPI.h"
#include "Mesh.h"
#include "Potentials.h"
#include "Rho.h"
#include "XConGrid.h"

#include <xc.h>

#include <vector>

template <class T>
class LDAonGridLibXC : public XConGrid
{
    Rho<T>& rho_;

    xc_func_type xfunc_;
    xc_func_type cfunc_;
    std::vector<double> exc_;
    std::vector<double> vxc_;

    Potentials& pot_;

public:
    LDAonGridLibXC(Rho<T>& rho, Potentials& pot) : rho_(rho), pot_(pot)
    {
        int func_id = XC_LDA_X;
        if (xc_func_init(&xfunc_, func_id, XC_UNPOLARIZED) != 0)
        {
            std::cerr << "Functional " << func_id << " not found" << std::endl;
        }
        func_id = XC_LDA_C_PZ_MOD;
        if (xc_func_init(&cfunc_, func_id, XC_UNPOLARIZED) != 0)
        {
            std::cerr << "Functional " << func_id << " not found" << std::endl;
        }
        exc_.resize(rho.rho_[0].size());
        vxc_.resize(rho.rho_[0].size());
    }

    ~LDAonGridLibXC() override
    {
        xc_func_end(&xfunc_);
        xc_func_end(&cfunc_);
    }

    void update() override;

    double getExc() const override; // in [Ha]
};

#endif

#endif
