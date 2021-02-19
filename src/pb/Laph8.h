// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PB_LAPH8_H
#define PB_LAPH8_H

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
    using memory_space_type = FDoperInterface::memory_space_type;

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
        lower_order_op_ = nullptr;
    }

    ~Laph8() override
    {
        if (lower_order_op_ != nullptr)
        {
            delete lower_order_op_;
            lower_order_op_ = nullptr;
        }
    }

    // construct a coarse grid operator
    Laph8 coarseOp(const Grid& mygrid)
    {
        Grid coarse_G(mygrid.coarse_grid(), 4);

        Laph8 A(coarse_G);

        return A;
    }

    Laph8 replicatedOp(const Grid& replicated_grid)
    {
        Laph8 replicated_A(replicated_grid);

        return replicated_A;
    }

    void setLowerOrderGrid() override
    {
        FDoper<T>::setFDLowerOrderGrid(Laph6<T>::minNumberGhosts());
    }

    Laph6<T>& getLowerOrderOp()
    {
        if (lower_order_op_ == nullptr)
        {
            FDoper<T>::setFDLowerOrderGrid(Laph6<T>::minNumberGhosts());
            lower_order_op_ = new Laph6<T>(Lap<T>::getLowerOrderGrid());
        }
        return *lower_order_op_;
    }

    static int minNumberGhosts() { return 4; }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B) override
    {
        // if(grid_.mype_env().mytask()==0)cout<<Lap<T>::name_<<endl;
        FDoper<T>::del2_8th(A, B);
        B.set_bc(A.bc(0), A.bc(1), A.bc(2));
    }

    void jacobi(GridFunc<T>&, const GridFunc<T>&, GridFunc<T>&) override;

    double jacobiFactor() const override { return invDiagEl_ / 1.5; }

    double diagEl(void) const override { return diagEl_; };
    double invDiagEl(void) const override { return invDiagEl_; };
};

} // namespace pb

#endif
