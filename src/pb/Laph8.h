// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Laph8.h,v 1.6 2009/07/16 23:36:28 jeanluc Exp $
#ifndef LAPH8_H
#define LAPH8_H

#include "Laph6.h"

// Laplacian operator O(h^8)
namespace pb
{
template <class T>
class Laph8 : public Lap<T>
{
    double diagEl_;
    double invDiagEl_;

    Laph6<T>* lower_order_op_;

public:
    Laph8(const Grid& mygrid) : Lap<T>(mygrid)
    {
        // cout<<" Create Laph8 operator\n";

        assert(fabs((Lap<T>::inv_h2(0) - Lap<T>::inv_h2(1)) / Lap<T>::inv_h2(0))
               < 0.05);
        assert(fabs((Lap<T>::inv_h2(1) - Lap<T>::inv_h2(2)) / Lap<T>::inv_h2(1))
               < 0.05);
        assert(fabs((Lap<T>::inv_h2(2) - Lap<T>::inv_h2(0)) / Lap<T>::inv_h2(2))
               < 0.05);

        diagEl_ = (1435. / 504.)
                  * (Lap<T>::inv_h2(0) + Lap<T>::inv_h2(1) + Lap<T>::inv_h2(2));
        invDiagEl_      = 1. / diagEl_;
        Lap<T>::name_   = "8th order";
        lower_order_op_ = NULL;
    }

    ~Laph8()
    {
        if (lower_order_op_ != NULL)
        {
            delete lower_order_op_;
            lower_order_op_ = NULL;
        }
    }

    // construct a coarse grid operator
    Laph8 coarseOp(const Grid& mygrid)
    {
        Grid coarse_G = mygrid.coarse_grid();

        Laph8 A(coarse_G);

        return A;
    }

    Laph8 replicatedOp(const Grid& replicated_grid)
    {
        Laph8 replicated_A(replicated_grid);

        return replicated_A;
    }

    void setLowerOrderGrid()
    {
        FDoper<T>::setFDLowerOrderGrid(Laph6<T>::minNumberGhosts());
    }

    Laph6<T>& getLowerOrderOp()
    {
        if (lower_order_op_ == NULL)
        {
            FDoper<T>::setFDLowerOrderGrid(Laph6<T>::minNumberGhosts());
            lower_order_op_ = new Laph6<T>(Lap<T>::getLowerOrderGrid());
        }
        return *lower_order_op_;
    }

    static int minNumberGhosts() { return 4; }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B)
    {
        // if(grid_.mype_env().mytask()==0)cout<<Lap<T>::name_<<endl;
        FDoper<T>::del2_8th(A, B);
        B.set_bc(A.bc(0), A.bc(1), A.bc(2));
    }
    void apply(GridFuncVector<T>& A, GridFuncVector<T>& B)
    {
        assert(A.size() == B.size());
        A.trade_boundaries();
        const int nfunc = (int)A.size();
        for (int k = 0; k < nfunc; k++)
        {
            FDoper<T>::del2_8th(A.func(k), B.func(k));
        }
    }

    void jacobi(GridFunc<T>&, const GridFunc<T>&, GridFunc<T>&);
    void jacobi(GridFuncVector<T>&, const GridFuncVector<T>&, GridFunc<T>&);
    void jacobi(
        GridFuncVector<T>&, const GridFuncVector<T>&, GridFuncVector<T>&);

    double diagEl(void) const { return diagEl_; };
    double invDiagEl(void) const { return invDiagEl_; };
};

} // namespace pb

#endif
