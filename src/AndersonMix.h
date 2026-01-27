// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ANDERSON_MIX
#define MGMOL_ANDERSON_MIX

#include "Mixing.h"
#include "Timer.h"

#include <cassert>
#include <iostream>
#include <vector>

template <class T>
class AndersonMix : public Mixing<T>
{
    const int m_;
    int mm_;
    double beta_; // mixing parameter

    std::vector<T*> xi_; // last mm_ trial solutions
    std::vector<T*> fi_; // last mm_ residuals
    std::vector<double> mat_;
    std::vector<double> rhs_;
    std::vector<double> theta_;

    static Timer update_tm_;

    virtual void postprocessUpdate() {};

    T& x_; // current trial solution

public:
    static Timer update_tm() { return update_tm_; }

    AndersonMix(const int m, const double beta, T& x);

    ~AndersonMix() override;

    // update trial solution based on residual
    // need work array for temporary storage
    void update(T& res, T& work, std::ostream& os, const bool verbose) override;
    void restart(void) override;
};

template <class T>
Timer AndersonMix<T>::update_tm_("AndersonMix<T>::update");

#endif
