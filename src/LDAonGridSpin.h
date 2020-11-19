// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_LDAONGRIDSPIN_H
#define MGMOL_LDAONGRIDSPIN_H

#include "LDAFunctional.h"
#include "MGmol_MPI.h"
#include "Mesh.h"
#include "Rho.h"
#include "XConGrid.h"

//#define USE_LIBXC

#ifdef USE_LIBXC
#include "Control.h"
#include <xc.h>
#endif

#include <vector>

class Potentials;

template <class T>
class LDAonGridSpin : public XConGrid
{
    int np_;
    int myspin_;

#ifdef USE_LIBXC
    xc_func_type xfunc_;
    xc_func_type cfunc_;
    std::vector<double> exc_;
    std::vector<double> vxc_;
#else
    LDAFunctional* lda_;
#endif
    Rho<T>& rho_;

    Potentials& pot_;

public:
    LDAonGridSpin(Rho<T>& rho, Potentials& pot)
        : np_(rho.rho_[0].size()), rho_(rho), pot_(pot)
    {
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        myspin_         = mmpi.myspin();
#ifdef USE_LIBXC
        int func_id = XC_LDA_X;
        if (xc_func_init(&xfunc_, func_id, XC_POLARIZED) != 0)
        {
            cerr << "Functional " << func_id << " not found" << endl;
        }
        func_id = XC_LDA_C_PZ_MOD;
        if (xc_func_init(&cfunc_, func_id, XC_POLARIZED) != 0)
        {
            cerr << "Functional " << func_id << " not found" << endl;
        }
        exc_.resize(np_ * 2);
        vxc_.resize(np_ * 2);
#else
        lda_ = new LDAFunctional(rho.rho_);
#endif
    }

    ~LDAonGridSpin() override
    {
#ifdef USE_LIBXC
        xc_func_end(&xfunc_);
        xc_func_end(&cfunc_);
#else
        delete lda_;
#endif
    }

    void update() override;

    double getExc() const override // in [Ha]
    {
        Mesh* mymesh           = Mesh::instance();
        const pb::Grid& mygrid = mymesh->grid();

#ifdef USE_LIBXC
        double sum = 0.;
        double exc = mygrid.vel() * MPdot(np_, &rho_.rho_[0][0], &exc_[0]);
        exc += ddot(&np_, &rho_.rho_[1][0], &ione, &exc_[0], &ione);
        //        double exc= mygrid.vel()*ddot(&np, &rho_.rho_[0][0], &ione,
        //        &exc_[0], &ione);

        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        int rc          = mmpi.allreduce(&exc, &sum, 1, MPI_SUM);
        if (rc != MPI_SUCCESS)
        {
            (*MPIdata::sout) << "MPI_Allreduce double sum failed!!!" << endl;
            Control& ct = *(Control::instance());
            ct.global_exit(2);
        }
        return sum;
#else
        return mygrid.vel() * lda_->computeRhoDotExc();
#endif
    }
};

#endif
