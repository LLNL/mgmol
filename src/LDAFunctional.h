// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef LDAFUNCTIONAL_H
#define LDAFUNCTIONAL_H

#include "Rho.h"
#include "XCFunctional.h"

#include <cassert>
#include <string>
#include <vector>

class LDAFunctional : public XCFunctional
{
    void xc_unpolarized(RHODTYPE rh, POTDTYPE& ee, POTDTYPE& vv);
    void xc_polarized(RHODTYPE rh, POTDTYPE& ee, POTDTYPE& vv);
    std::vector<POTDTYPE> exc_;
    std::vector<std::vector<POTDTYPE>> vxc_;

    LDAFunctional();

public:
    LDAFunctional(std::vector<std::vector<RHODTYPE>>& rhoe) : XCFunctional(rhoe)
    {
        assert(np_ > 0);

        exc_.resize(np_);
        vxc_.resize(nspin_);
        for (int i = 0; i < nspin_; i++)
        {
            vxc_[i].resize(np_);
        }

        pexc_ = &exc_[0];
        if (nspin_ == 1)
        {
            pvxc1_ = &vxc_[0][0];
        }
        else
        {
            pvxc1_up_ = &vxc_[0][0];
            pvxc1_dn_ = &vxc_[1][0];
        }
    };

    bool isGGA() const override { return false; };
    std::string name() const override { return "LDA"; };
    void computeXC(void) override;

    double computeRhoDotExc() const;
};
#endif
