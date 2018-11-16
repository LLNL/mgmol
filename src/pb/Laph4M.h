// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Laph4M.h,v 1.12 2009/07/17 00:07:45 jeanluc Exp $
#ifndef PB_LAPH4M_H
#define PB_LAPH4M_H

#include "Lap.h"
namespace pb
{
template <class T>
class Laph4M : public Lap<T>
{
    double diagEl_;
    double invDiagEl_;

public:
    Laph4M(const Grid& mygrid) : Lap<T>(mygrid)
    {
        // cout<<" Create Laph4M operator\n";

        diagEl_ = (4. / 3.)
                  * (Lap<T>::inv_h2(0) + Lap<T>::inv_h2(1) + Lap<T>::inv_h2(2));
        invDiagEl_    = 1. / diagEl_;
        Lap<T>::name_ = "Mehrstellen 4th order";
    }

    // construct a coarse grid operator
    Laph4M coarseOp(const Grid& mygrid)
    {
        Grid coarse_G = mygrid.coarse_grid();

        Laph4M A(coarse_G);

        return A;
    }

    Laph4M replicatedOp(const Grid& replicated_grid)
    {
        Laph4M replicated_A(replicated_grid);

        return replicated_A;
    }

    void setLowerOrderGrid()
    {
        FDoper<T>::setFDLowerOrderGrid(Laph4M::minNumberGhosts());
    }

    Laph4M& getLowerOrderOp() { return *this; }

    static int minNumberGhosts() { return 1; }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B)
    {
        FDoper<T>::del2_4th_Mehr(A, B);
        B.set_bc(A.bc(0), A.bc(1), A.bc(2));
    }
    void apply(GridFuncVector<T>& A, GridFuncVector<T>& B)
    {
        assert(A.size() == B.size());
        A.trade_boundaries();
        const int nfunc = (int)A.size();
        for (int k = 0; k < nfunc; k++)
        {
            Lap<T>::del2_4th_Mehr(A.func(k), B.func(k));
        }
    }

    void rhs(GridFunc<T>& A, GridFunc<T>& B) const
    {
        this->rhs_4th_Mehr1(A, B);
        B.set_bc(A.bc(0), A.bc(1), A.bc(2));
    }
    void rhs(GridFunc<T>& A, T* const B) const { this->rhs_4th_Mehr1(A, B); }

    void jacobi(GridFunc<T>&, const GridFunc<T>&, GridFunc<T>&);
    void jacobi(GridFuncVector<T>&, const GridFuncVector<T>&, GridFunc<T>&);
    void jacobi(
        GridFuncVector<T>&, const GridFuncVector<T>&, GridFuncVector<T>&);

    double diagEl(void) const { return diagEl_; };
    double invDiagEl(void) const { return invDiagEl_; };
};

} // namespace pb

#endif
