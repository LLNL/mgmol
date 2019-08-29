// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "LDAonGrid.h"
#include "PBEonGrid.h"
#include "PBEonGridSpin.h"
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
            return new LDAonGrid<T>(rho, pot);
        }
        else if (xctype == 2)
        {
            if (nspin > 1)
                return new PBEonGridSpin<T>(rho, pot);
            else
                return new PBEonGrid<T>(rho, pot);
        }
        else
        {
            std::cerr << "Invalid XC option" << std::endl;
            return nullptr;
        }
    }
};
