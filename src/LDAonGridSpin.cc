// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "LDAonGridSpin.h"

#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Potentials.h"

template <class T>
void LDAonGridSpin<T>::update()
{
    get_xc_tm_.start();

    int iterative_index = rho_.getIterativeIndex();

    lda_->computeXC();
    if (myspin_ == 0)
    {
        pot_.setVxc(lda_->pvxc1_up_, iterative_index);
    }
    else
    {
        pot_.setVxc(lda_->pvxc1_dn_, iterative_index);
    }

    get_xc_tm_.stop();
}

template <class T>
double LDAonGridSpin<T>::getExc() const // in [Ha]
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    return mygrid.vel() * lda_->computeRhoDotExc();
}

template class LDAonGridSpin<LocGridOrbitals<ORBDTYPE>>;
template class LDAonGridSpin<ExtendedGridOrbitals<ORBDTYPE>>;
