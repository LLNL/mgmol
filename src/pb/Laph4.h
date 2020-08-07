// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PB_LAPH4_H
#define PB_LAPH4_H

#include "Laph2.h"

// Laplacian operator O(h^4)
namespace pb
{
template <class T>
class Laph4 : public Lap<T>
{
    double diagEl_;
    double invDiagEl_;

    Laph2<T>* lower_order_op_;

public:
    Laph4(const Grid& mygrid) : Lap<T>(mygrid)
    {
        // cout<<" Create Laph4 operator\n";

        //-1/12 16/12 -30/12 16/12 -1/12
        diagEl_ = 2.5
                  * (Lap<T>::inv_h2(0) + Lap<T>::inv_h2(1)
                        + Lap<T>::inv_h2(2)); // 2.5 = 30./12.
        invDiagEl_    = 1. / diagEl_;
        Lap<T>::name_ = "Classical 4th order";

        lower_order_op_ = nullptr;
    }

    ~Laph4() override
    {
        if (lower_order_op_ != nullptr)
        {
            delete lower_order_op_;
            lower_order_op_ = nullptr;
        }
    }

    // construct a coarse grid operator
    Laph4 coarseOp(const Grid& mygrid)
    {
        Grid coarse_G(mygrid.coarse_grid(), 2);

        Laph4 A(coarse_G);

        return A;
    }

    Laph4 replicatedOp(const Grid& replicated_grid)
    {
        Laph4 replicated_A(replicated_grid);

        return replicated_A;
    }

    void setLowerOrderGrid() override
    {
        this->setFDLowerOrderGrid(Laph2<T>::minNumberGhosts());
    }

    Laph2<T>& getLowerOrderOp()
    {
        if (lower_order_op_ == nullptr)
        {
            this->setFDLowerOrderGrid(Laph2<T>::minNumberGhosts());
            lower_order_op_ = new Laph2<T>(Lap<T>::getLowerOrderGrid());
        }
        return *lower_order_op_;
    }

    static short minNumberGhosts() { return 2; }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B) override
    {
        if (!A.updated_boundaries()) A.trade_boundaries();
        this->del2_4th(A.grid(), A.uu(), B.uu(), 1);
        B.set_updated_boundaries(0);
        B.set_bc(A.bc(0), A.bc(1), A.bc(2));
    }
    void applyWithPot(GridFunc<T>& A, const double* pot, T* B) override
    {
        this->del2_4th_withPot(A, pot, B);
    }
    void apply(GridFuncVector<T>& A, GridFuncVector<T>& B) override
    {
        assert(A.size() == B.size());
        A.trade_boundaries();
        const int nfunc = (int)A.size();
        for (int k = 0; k < nfunc; k++)
        {
            this->del2_4th(A.getGridFunc(k), B.getGridFunc(k));
        }
    }
    void apply(Grid& Agrid, T* A, T* B, const size_t nfunc)
    {
        this->del2_4th(Agrid, A, B, nfunc);
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
