// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "LDAonGrid.h"
#include "LDAonGridLibXC.h"
#include "LDAonGridSpin.h"
#include "LDAonGridSpinLibXC.h"
#include "PBEonGrid.h"
#include "PBEonGridLibXC.h"
#include "PBEonGridSpin.h"
#include "PBEonGridSpinLibXC.h"
#include "Potentials.h"

#include <iostream>

template <class T>
class XCfunctionalFactory
{
public:
    static XConGrid* create(
        const int xctype, const int nspin, Rho<T>& rho, Potentials& pot)
    {
        if (xctype == 0)
        {
            if (nspin > 1)
#ifdef USE_LIBXC
                return new LDAonGridSpinLibXC<T>(rho, pot);
#else
                return new LDAonGridSpin<T>(rho, pot);
#endif
            else
#ifdef USE_LIBXC
                return new LDAonGridLibXC<T>(rho, pot);
#else
                return new LDAonGrid<T>(rho, pot);
#endif
        }
        else if (xctype == 2)
        {
            if (nspin > 1)
#ifdef USE_LIBXC
                return new PBEonGridSpinLibXC<T>(rho, pot);
#else
                return new PBEonGridSpin<T>(rho, pot);
#endif
            else
#ifdef USE_LIBXC
                return new PBEonGridLibXC<T>(rho, pot);
#else
                return new PBEonGrid<T>(rho, pot);
#endif
        }
        else
        {
            std::cerr << "Invalid XC option" << std::endl;
            return nullptr;
        }
    }
};
