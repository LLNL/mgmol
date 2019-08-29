// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Laph2.h,v 1.12 2010/01/28 22:56:47 jeanluc Exp $
#ifndef PB_LAPH2_H
#define PB_LAPH2_H

#include "Lap.h"

namespace pb
{
template <class T>
class Laph2 : public Lap<T>
{
    double diagEl_;
    double invDiagEl_;

    Laph2* lower_order_op_;

public:
    Laph2(const Grid& mygrid) : Lap<T>(mygrid)
    {
        // cout<<" Create Laph2 operator\n";

        diagEl_
            = 2. * (Lap<T>::inv_h2(0) + Lap<T>::inv_h2(1) + Lap<T>::inv_h2(2));
        invDiagEl_    = 1. / diagEl_;
        Lap<T>::name_ = "2nd order";

        lower_order_op_ = nullptr;
    }

    ~Laph2() override
    {
        if (lower_order_op_ != nullptr)
        {
            delete lower_order_op_;
            lower_order_op_ = nullptr;
        }
    }

    // construct a coarse grid operator
    Laph2 coarseOp(const Grid& mygrid)
    {
        Grid coarse_G = mygrid.coarse_grid();

        Laph2 A(coarse_G);

        return A;
    }

    Laph2 replicatedOp(const Grid& replicated_grid) const
    {
        Laph2 replicated_A(replicated_grid);

        return replicated_A;
    }

    void setLowerOrderGrid() override
    {
        this->setFDLowerOrderGrid(minNumberGhosts());
    }

    Laph2& getLowerOrderOp()
    {
        if (lower_order_op_ == nullptr)
        {
            this->setFDLowerOrderGrid(Laph2::minNumberGhosts());
            lower_order_op_ = new Laph2(Lap<T>::getLowerOrderGrid());
        }
        return *lower_order_op_;
    }

    static short minNumberGhosts() { return 1; }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B) override
    {
        // if(grid_.mype_env().mytask()==0)cout<<Lap<T>::name_<<endl;
        this->del2_2nd(A, B);
        B.set_bc(A.bc(0), A.bc(1), A.bc(2));
    }
    void apply(GridFuncVector<T>& A, GridFuncVector<T>& B) override
    {
        assert(A.size() == B.size());
        A.trade_boundaries();
        const int nfunc = (int)A.size();
        for (int k = 0; k < nfunc; k++)
        {
            this->del2_2nd(A.func(k), B.func(k));
        }
    }

    void jacobi(GridFunc<T>&, const GridFunc<T>&, GridFunc<T>&) override;
    void jacobi(
        GridFuncVector<T>&, const GridFuncVector<T>&, GridFunc<T>&) override;
    void jacobi(GridFuncVector<T>&, const GridFuncVector<T>&,
        GridFuncVector<T>&) override;

    double diagEl(void) const override { return diagEl_; };
    double invDiagEl(void) const override { return invDiagEl_; };
};

} // namespace pb

#endif
