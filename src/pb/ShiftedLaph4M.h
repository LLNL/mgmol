// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: ShiftedLaph4M.h,v 1.12 2009/07/17 00:07:45 jeanluc Exp $
#ifndef ShiftedLaph4M_H
#define ShiftedLaph4M_H

#include "Lap.h"

namespace pb
{
template <class T>
class ShiftedLaph4M : public Lap<T>
{

private:
    double diagEl_;
    double invDiagEl_;
    double lambda2_;

public:
    ShiftedLaph4M(const Grid& mygrid, const double shift) : Lap<T>(mygrid)
    {
        // cout<<" Create ShiftedLaph4M operator\n";
        lambda2_ = shift;
        if (lambda2_ < 0.)
        {
            std::cout
                << "Warning! ShiftedLaph4M constructor: shift should be >0"
                << std::endl;
        }
        diagEl_
            = (4. / 3.)
                  * (Lap<T>::inv_h2(0) + Lap<T>::inv_h2(1) + Lap<T>::inv_h2(2))
              + lambda2_;
        invDiagEl_    = 1. / diagEl_;
        Lap<T>::name_ = "Shifted Mehrstellen 4th order";
    }

    // construct a coarse grid operator
    ShiftedLaph4M coarseOp(const Grid& mygrid)
    {
        Grid coarse_G = mygrid.coarse_grid();

        ShiftedLaph4M A(coarse_G, lambda2_);

        return A;
    }

    ShiftedLaph4M replicatedOp(const Grid& replicated_grid)
    {
        ShiftedLaph4M replicated_A(replicated_grid, lambda2_);

        return replicated_A;
    }

    void setLowerOrderGrid()
    {
        this->setFDLowerOrderGrid(ShiftedLaph4M::minNumberGhosts());
    }

    ShiftedLaph4M& getLowerOrderOp() { return *this; }

    static int minNumberGhosts() { return 1; }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B);
    void apply(GridFuncVector<T>& A, GridFuncVector<T>& B);

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

    double lambda2() const { return lambda2_; }

    double diagEl(void) const { return diagEl_; };
    double invDiagEl(void) const { return invDiagEl_; };
};

} // namespace pb

#endif
