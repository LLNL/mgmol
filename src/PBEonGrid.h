// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef PBEONGRID_H
#define PBEONGRID_H

#include "PBEFunctional.h"
#include "Rho.h"
#include "XConGrid.h"
#include "Mesh.h"

//#define USE_LIBXC

#ifdef USE_LIBXC
#include <xc.h>
#include "MGmol_MPI.h"
#include "Control.h"
#endif

#include <vector>

class Potentials;

class PBEonGrid : public XConGrid
{
    int np_;
#ifdef USE_LIBXC
    xc_func_type xfunc_;
    xc_func_type cfunc_;
    std::vector<double> exc_;
    std::vector<double> vxc_;
    std::vector<double> vsigma_;
#else
    PBEFunctional* pbe_;
#endif
    Rho& rho_;
    
    Potentials& pot_;

public:
    PBEonGrid(Rho& rho, Potentials& pot):
        rho_(rho),pot_(pot)
    {
        np_=rho.rho_[0].size();
#ifdef USE_LIBXC
        int func_id = XC_GGA_X_PBE;
        if(xc_func_init(&xfunc_, func_id, XC_UNPOLARIZED) != 0){
            cerr<<"Functional "<<func_id<<" not found"<<endl;
        }
        func_id = XC_GGA_C_PBE;
        if(xc_func_init(&cfunc_, func_id, XC_UNPOLARIZED) != 0){
            cerr<<"Functional "<<func_id<<" not found"<<endl;
        }
        exc_.resize(np_);
        vxc_.resize(np_);
        vsigma_.resize(np_);
#else
        pbe_=new PBEFunctional(rho.rho_);
#endif
    }
    
    ~PBEonGrid()
    {
#ifdef USE_LIBXC
        xc_func_end (&xfunc_);
        xc_func_end (&cfunc_);
#else
        delete pbe_;
#endif
    }
    
    void update();

    double getExc()const;
};

#endif
