// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: PB.cc,v 1.2 2008/12/10 00:58:56 jeanluc Exp $
#include "PB.h"

namespace pb
{
template <class T>
double PB<T>::energyES(GridFunc<T>& v, GridFunc<T>& rho)
{
    double g = 0.5 * dot(v, rho);

    if (FDoper<T>::grid_.mype_env().mytask() == 0)
        std::cout << " ES Energy = " << g << std::endl;

    return g;
}
template class PB<double>;
template class PB<float>;
} // namespace pb
