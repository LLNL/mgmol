// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PB_FDOPER_H
#define PB_FDOPER_H

#include "FDoperInterface.h"

#include "Grid.h"
#include "GridFunc.h"
#include "memory_space.h"

namespace pb
{

template <class T>
class FDoper : public FDoperInterface
{
    Grid* lower_order_grid_;

protected:
    Grid grid_;

    double inv_h_[3];
    double inv_h2_[3];
    int dim_[3];
    int incx_;
    int incy_;

    void del1_4th(GridFunc<T>&, GridFunc<T>&, const short) const;
    void del2_4th(GridFunc<T>&, GridFunc<T>&) const;
    void del2_4th_withPot(GridFunc<T>&, const double* const pot, T*) const;
    void del1_2nd(GridFunc<T>&, GridFunc<T>&, const short) const;
    void del2_2nd(GridFunc<T>&, GridFunc<T>&) const;
    void del1_6th(GridFunc<T>&, GridFunc<T>&, const short) const;
    void del2_6th(GridFunc<T>&, GridFunc<T>&) const;
    void del1_8th(GridFunc<T>&, GridFunc<T>&, const short) const;
    void del2_8th(GridFunc<T>&, GridFunc<T>&) const;

    // Mehrstellenverfahren operators
    void del2_4th_Mehr(GridFunc<T>&, GridFunc<T>&) const;
    void rhs_4th_Mehr1(GridFunc<T>&, GridFunc<T>&) const;
    void rhs_4th_Mehr1(GridFunc<T>&, T* const) const;
    void rhs_4th_Mehr2(GridFunc<T>&, GridFunc<T>&) const;
    void rhs_4th_Mehr2(GridFunc<T>&, T* const) const;
    double Mgm(GridFunc<T>&, const GridFunc<T>&, const int, const int) const;

public:
    FDoper(const Grid&);
    FDoper(const FDoper& oper, const Grid&);

    int dim(const int i) const { return dim_[i]; }
    short ghosts() const { return grid_.ghost_pt(); };
    double inv_h(const short i) const { return inv_h_[i]; }
    double inv_h2(const short i) const { return inv_h2_[i]; }
    const Grid& grid() { return grid_; }

    void smooth(GridFunc<T>&, GridFunc<T>&, const double);

    virtual void transform(GridFunc<T>&) const {};
    virtual void inv_transform(GridFunc<T>&) const {};
    virtual void rhs(GridFunc<T>& A, GridFunc<T>& B) const
    {
        rhs_tm_.start();
        B = A;
        rhs_tm_.stop();
    };
    virtual void rhs(GridFunc<T>& A, T* B) const
    {
        rhs_tm_.start();
        A.init_vect(B, 'd');
        rhs_tm_.stop();
    }

    const Grid& getLowerOrderGrid() const
    {
        assert(lower_order_grid_);
        return *lower_order_grid_;
    }

    void setFDLowerOrderGrid(const short nghosts)
    {
        assert(nghosts >= 1);
        if (lower_order_grid_ == nullptr)
            lower_order_grid_ = new Grid(grid_, nghosts);
    }

    virtual void apply(GridFunc<T>&, GridFunc<T>&) = 0;

    GridFunc<T> operator*(GridFunc<T>& A)
    {
        GridFunc<T> work(A.grid(), A.bc(0), A.bc(1), A.bc(2));

        apply(A, work);

        return work;
    }

    ~FDoper() override
    {
        if (lower_order_grid_ != nullptr)
        {
            delete lower_order_grid_;
            lower_order_grid_ = nullptr;
        }
    }
};

} // namespace pb

#endif
