// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Laph6.h,v 1.11 2010/01/28 22:56:47 jeanluc Exp $
#ifndef LAPH6_H
#define LAPH6_H

#include "Laph4.h"

// Laplacian operator O(h^6)
namespace pb
{
template <class T>
class Laph6 : public Lap<T>
{
    double diagEl_;
    double invDiagEl_;

    Laph4<T>* lower_order_op_;

public:
    Laph6(const Grid& mygrid) : Lap<T>(mygrid)
    {
        // cout<<" Create Laph6 operator\n";

        assert(fabs((Lap<T>::inv_h2(0) - Lap<T>::inv_h2(1)) / Lap<T>::inv_h2(0))
               < 0.05);
        assert(fabs((Lap<T>::inv_h2(1) - Lap<T>::inv_h2(2)) / Lap<T>::inv_h2(1))
               < 0.05);
        assert(fabs((Lap<T>::inv_h2(2) - Lap<T>::inv_h2(0)) / Lap<T>::inv_h2(2))
               < 0.05);

        diagEl_ = (49. / 18.)
                  * (Lap<T>::inv_h2(0) + Lap<T>::inv_h2(1) + Lap<T>::inv_h2(2));
        invDiagEl_    = 1. / diagEl_;
        Lap<T>::name_ = "6th order";

        lower_order_op_ = nullptr;
    }

    ~Laph6() override
    {
        if (lower_order_op_ != nullptr)
        {
            delete lower_order_op_;
            lower_order_op_ = nullptr;
        }
    }

    // construct a coarse grid operator
    Laph6 coarseOp(const Grid& mygrid)
    {
        Grid coarse_G = mygrid.coarse_grid();

        Laph6 A(coarse_G);

        return A;
    }

    Laph6 replicatedOp(const Grid& replicated_grid)
    {
        Laph6 replicated_A(replicated_grid);

        return replicated_A;
    }

    void setLowerOrderGrid() override
    {
        FDoper<T>::setFDLowerOrderGrid(Laph4<T>::minNumberGhosts());
    }

    Laph4<T>& getLowerOrderOp()
    {
        if (lower_order_op_ == nullptr)
        {
            FDoper<T>::setFDLowerOrderGrid(Laph4<T>::minNumberGhosts());
            lower_order_op_ = new Laph4<T>(Lap<T>::getLowerOrderGrid());
        }
        return *lower_order_op_;
    }

    static short minNumberGhosts() { return 3; }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B) override
    {
        // if(grid_.mype_env().mytask()==0)cout<<Lap<T>::name_<<endl;
        FDoper<T>::del2_6th(A, B);
        B.set_bc(A.bc(0), A.bc(1), A.bc(2));
    }
    void apply(GridFuncVector<T>& A, GridFuncVector<T>& B) override
    {
        assert(A.size() == B.size());
        A.trade_boundaries();
        const int nfunc = (int)A.size();
        for (int k = 0; k < nfunc; k++)
        {
            FDoper<T>::del2_6th(A.func(k), B.func(k));
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
