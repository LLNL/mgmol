// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

////////////////////////////////////////////////////////////////////////////////
//
// Abstract base class for density functionals
// Input variables are: rho, rho_up, rho_dn, grad_rho, grad_rho_up, grad_rho_dn
//
// Output quantities:
// The exchange-correlation energy is expressed as
//
// Exc = int{ rho[i] * exc[i] + rho_up[i] * exc_up[i] + rho_dn[i] * exc_dn[i] }
//
// It is assumed that exchange correlation potentials can be
// written as:
//
// vxc_up = vxc1 + vxc1_up +
//          div ( vxc2_upup grad_rho_up ) + div ( vxc2_updn grad_rho_dn )
//
// vxc_dn = vxc1 + vxc1_dn +
//          div ( vxc2_dndn grad_rho_dn ) + div ( vxc2_dnup grad_rho_up )
//
// Not all input quantities are needed, and not all output quantities are
// computed by certain functionals. Example:
//
// LDAFunctional:
//   without spin: Input: rho
//                 Output: exc, vxc1
//   with spin:    Input: rho_up, rho_dn
//                 Output: vxc1_up, vxc1_dn
//
// PBEFunctional:
//   without spin: Input:  rho, grad_rho
//                 Output: exc, vxc1, vxc2
//   with spin:    Input:  rho_up, rho_dn, grad_rho_up, grad_rho_dn,
//                 Output: exc_up, exc_dn, vxc1_up, vxc1_dn,
//                         vxc2_upup, vxc2_dndn, vxc2_updn, vxc2_dnup

#ifndef XCFUNCTIONAL_H
#define XCFUNCTIONAL_H

#include "MGmol_blas1.h"

#include <cassert>
#include <string>
#include <vector>

#include "Rho.h"

class XCFunctional
{
protected:
    int np_;
    int nspin_;

    RHODTYPE* prho_;
    RHODTYPE* prho_up_;
    RHODTYPE* prho_dn_;

    POTDTYPE* pexc_;
    POTDTYPE* pexc_up_;
    POTDTYPE* pexc_dn_;

public:
    POTDTYPE *pvxc1_, *pvxc1_up_, *pvxc1_dn_;
    POTDTYPE *pvxc2_, *pvxc2_upup_, *pvxc2_dndn_, *pvxc2_updn_, *pvxc2_dnup_;

    virtual bool isGGA(void) const       = 0;
    virtual std::string name(void) const = 0;
    int np(void) const { return np_; };
    int nspin(void) const { return nspin_; };

    XCFunctional(std::vector<std::vector<RHODTYPE>>& rhoe)
        : np_(rhoe[0].size()), nspin_(rhoe.size())
    {
        if (nspin_ > 1) assert(rhoe[0].size() == rhoe[1].size());

        pexc_ = pexc_up_ = pexc_dn_ = nullptr;
        pvxc1_ = pvxc1_up_ = pvxc1_dn_ = nullptr;
        pvxc2_ = pvxc2_upup_ = pvxc2_dndn_ = pvxc2_updn_ = pvxc2_dnup_
            = nullptr;

        if (nspin_ == 1)
        {
            prho_ = &rhoe[0][0];

            prho_up_ = nullptr;
            prho_dn_ = nullptr;
        }
        else
        {
            prho_up_ = &rhoe[0][0];
            prho_dn_ = &rhoe[1][0];

            prho_ = nullptr;
        }
    }

    // virtual destructor needed to ensure proper deallocation
    virtual ~XCFunctional() { }

    virtual void computeXC(void) = 0;
};
#endif
