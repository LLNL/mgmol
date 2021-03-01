// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_LDAONGRID_H
#define MGMOL_LDAONGRID_H

#include "LDAFunctional.h"
#include "Mesh.h"
#include "Rho.h"
#include "XConGrid.h"

#include <vector>

class Potentials;

template <class T>
class LDAonGrid : public XConGrid
{
    Rho<T>& rho_;

    LDAFunctional* lda_;

    Potentials& pot_;

public:
    LDAonGrid(Rho<T>& rho, Potentials& pot) : rho_(rho), pot_(pot)
    {
        lda_ = new LDAFunctional(rho.rho_);
    }

    ~LDAonGrid() override { delete lda_; }

    void update() override;

    double getExc() const override // in [Ha]
    {
        Mesh* mymesh           = Mesh::instance();
        const pb::Grid& mygrid = mymesh->grid();

        return mygrid.vel() * lda_->computeRhoDotExc();
    }
};

#endif
