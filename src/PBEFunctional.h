// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PBEFUNCTIONAL_H
#define PBEFUNCTIONAL_H

#include "Rho.h"
#include "XCFunctional.h"

#include <string>
#include <vector>

class PBEFunctional : public XCFunctional
{
    PBEFunctional();

    std::vector<POTDTYPE> exc_;
    std::vector<POTDTYPE> exc_up_, exc_dn_;
    std::vector<POTDTYPE> vxc1_, vxc1_up_, vxc1_dn_, vxc2_, vxc2_upup_,
        vxc2_updn_, vxc2_dnup_, vxc2_dndn_;
    std::vector<RHODTYPE> grad_rho_[3], grad_rho_up_[3], grad_rho_dn_[3];

    void gcor2(double a, double a1, double b1, double b2, double b3, double b4,
        double rtrs, double* gg, double* ggrs);

    void excpbe(const RHODTYPE rho, const RHODTYPE grad, POTDTYPE* exc,
        POTDTYPE* vxc1, POTDTYPE* vxc2);

    void excpbe_sp(RHODTYPE rho_up, RHODTYPE rho_dn, RHODTYPE grad_up,
        RHODTYPE grad_dn, RHODTYPE grad, POTDTYPE* exc_up, POTDTYPE* exc_dn,
        POTDTYPE* vxc1_up, POTDTYPE* vxc1_dn, POTDTYPE* vxc2_upup,
        POTDTYPE* vxc2_dndn, POTDTYPE* vxc2_updn, POTDTYPE* vxc2_dnup);

    RHODTYPE* pgrad_rho_[3];
    RHODTYPE* pgrad_rho_up_[3];
    RHODTYPE* pgrad_rho_dn_[3];

public:
    PBEFunctional(std::vector<std::vector<RHODTYPE>>& rhoe);

    bool isGGA() const override { return true; };
    std::string name() const override { return "PBE"; };
    void computeXC(void) override;
    void setGradRho(const int i, const RHODTYPE* drho)
    {
        memcpy(pgrad_rho_[i], drho, np_ * sizeof(RHODTYPE));
    }
    void setGradRhoUp(const int i, const RHODTYPE* drho)
    {
        memcpy(pgrad_rho_up_[i], drho, np_ * sizeof(RHODTYPE));
    }
    void setGradRhoDn(const int i, const RHODTYPE* drho)
    {
        memcpy(pgrad_rho_dn_[i], drho, np_ * sizeof(RHODTYPE));
    }
    double computeRhoDotExc() const;
};
#endif
