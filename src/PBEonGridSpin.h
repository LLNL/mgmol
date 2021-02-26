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

#ifdef USE_LIBXC
#include <xc.h>
#else
#include "PBEFunctional.h"
#endif

#include <vector>

class Potentials;

template <class T>
class PBEonGridSpin : public XConGrid
{
    int np_;
    int myspin_;
    std::vector<double> vxc_;
#ifdef USE_LIBXC
    xc_func_type xfunc_;
    xc_func_type cfunc_;
    std::vector<double> exc_;
    std::vector<double> vsigma_;
#else
    PBEFunctional* pbe_;
#endif
    Rho<T>& rho_;

    Potentials& pot_;

public:
    PBEonGridSpin(Rho<T>& rho, Potentials& pot);

    ~PBEonGridSpin() override
    {
#ifdef USE_LIBXC
        xc_func_end(&xfunc_);
        xc_func_end(&cfunc_);
#else
        delete pbe_;
#endif
    }

    void update() override;

    double getExc() const override;
};

#endif
