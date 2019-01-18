// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ANDERSON_MIX
#define MGMOL_ANDERSON_MIX

#include "Mixing.h"
#include "Timer.h"

#include <cassert>
#include <vector>

template <class T>
class AndersonMix : public Mixing<T>
{
    const int m_;
    int mm_;
    double beta_; // mixing parameter
    T& x_; // current trial solution

    std::vector<T*> xi_; // last mm_ trial solutions
    std::vector<T*> fi_; // last mm_ residuals
    std::vector<double> mat_;
    std::vector<double> rhs_;
    std::vector<double> theta_;

    bool ortho_flag_;

    static Timer update_tm_;

public:
    static Timer update_tm() { return update_tm_; }

    AndersonMix(
        const int m, const double beta, T& x, const bool ortho_flag = false);

    ~AndersonMix();

    // update trial solution based on residual
    // need work array for temporary storage
    void update(T& res, T& work);
    void restart(void);
};

template <class T>
Timer AndersonMix<T>::update_tm_("AndersonMix<T>::update");

#endif
