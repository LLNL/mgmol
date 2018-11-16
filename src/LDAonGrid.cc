// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "LDAonGrid.h"
#include "Potentials.h"

void LDAonGrid::update()
{
    get_xc_tm_.start();

    int iterative_index = rho_.getIterativeIndex();

#ifdef USE_LIBXC
    xc_lda_exc_vxc(&xfunc_, exc_.size(), &rho_.rho_[0][0], &exc_[0], &vxc_[0]);
    vector<double> vtmp(vxc_.size());
    vector<double> etmp(vxc_.size());
    xc_lda_exc_vxc(&cfunc_, exc_.size(), &rho_.rho_[0][0], &etmp[0], &vtmp[0]);

    int ione   = 1;
    double one = 1.;
    int np     = vxc_.size();
    daxpy(&np, &one, &vtmp[0], &ione, &vxc_[0], &ione);
    daxpy(&np, &one, &etmp[0], &ione, &exc_[0], &ione);

    pot_.setVxc(&vxc_[0], iterative_index);
#else
    lda_->computeXC();
    pot_.setVxc(lda_->pvxc1_, iterative_index);
#endif

    get_xc_tm_.stop();
}
