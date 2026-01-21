// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "LDAonGrid.h"

#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Potentials.h"

template <class T>
void LDAonGrid<T>::update()
{
    get_xc_tm_.start();

    int iterative_index = rho_.getIterativeIndex();

    lda_->computeXC();
    pot_.setVxc(lda_->pvxc1_, iterative_index);

    get_xc_tm_.stop();
}

template class LDAonGrid<LocGridOrbitals<ORBDTYPE>>;
template class LDAonGrid<ExtendedGridOrbitals<ORBDTYPE>>;
