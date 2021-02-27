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
#include "memory_space.h"
#include "mputils.h"

//#define USE_LIBXC

#ifdef USE_LIBXC
#include "Control.h"
#include "MGmol_MPI.h"
#include <xc.h>
#endif

#include <vector>

class Potentials;

template <class T>
class LDAonGrid : public XConGrid
{
    Rho<T>& rho_;

#ifdef USE_LIBXC
    xc_func_type xfunc_;
    xc_func_type cfunc_;
    std::vector<double> exc_;
    std::vector<double> vxc_;
#else
    LDAFunctional* lda_;
#endif

    Potentials& pot_;

public:
    LDAonGrid(Rho<T>& rho, Potentials& pot) : rho_(rho), pot_(pot)
    {
#ifdef USE_LIBXC
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
#else
        lda_ = new LDAFunctional(rho.rho_);
#endif
    }

    ~LDAonGrid() override
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
        int np = exc_.size();
        //        int ione=1;
        double exc = mygrid.vel()
                     * LinearAlgebraUtils<MemorySpace::Host>::MPdot(
                           np, &rho_.rho_[0][0], &exc_[0]);

        double sum      = 0.;
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        int rc          = mmpi.allreduce(&exc, &sum, 1, MPI_SUM);
        if (rc != MPI_SUCCESS)
        {
            (*MPIdata::sout)
                << "MPI_Allreduce double sum failed!!!" << std::endl;
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
