// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: PBh6.h,v 1.11 2010/01/28 22:56:47 jeanluc Exp $
#ifndef PB_H6_H
#define PB_H6_H

#include "PB.h"
#include "PBh4.h"

// operator -epsilon div*grad -(grad epsilon)* grad
namespace pb
{
template <class T>
class PBh6 : public PB<T>
{

private:
    void pb_6th(GridFunc<T>& A, GridFunc<T>& B);

    PBh4<T>* lower_order_op_;

public:
    // constructor
    PBh6(const Grid& mygrid, DielFunc<T>& myepsilon) : PB<T>(mygrid, myepsilon)
    {
        PB<T>::initialized_ = true;
        lower_order_op_     = nullptr;
    };
    PBh6(const Grid& mygrid, const double e0, const double rho0,
        const double drho0)
        : PB<T>(mygrid, e0, rho0, drho0)
    {
        PB<T>::initialized_ = false;
        lower_order_op_     = nullptr;
    };

    ~PBh6()
    {
        if (lower_order_op_ != nullptr)
        {
            delete lower_order_op_;
            lower_order_op_ = nullptr;
        }
    }

    // construct a coarse grid operator
    PBh6 coarseOp(const Grid& mygrid)
    {
        Grid coarse_G = mygrid.coarse_grid();
        DielFunc<T> ecoarse(coarse_G, PB<T>::epsilon_.epsilon_max());
        PB<T>::epsilon_.restrict3D(ecoarse);

        PBh6 A(coarse_G, ecoarse);

        return A;
    }

    PBh6 replicatedOp(const Grid&);

    void setLowerOrderGrid()
    {
        FDoper<T>::setFDLowerOrderGrid(PBh4<T>::minNumberGhosts());
    }

    PBh4<T>& getLowerOrderOp()
    {
        if (lower_order_op_ == nullptr)
        {
            FDoper<T>::setFDLowerOrderGrid(PBh4<T>::minNumberGhosts());
            lower_order_op_
                = new PBh4<T>(PB<T>::getLowerOrderGrid(), PB<T>::epsilon_);
        }
        return *lower_order_op_;
    }

    static short minNumberGhosts() { return 3; }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B)
    {
        pb_6th(A, B);
        B.set_bc(A.bc(0), A.bc(1), A.bc(2));
    }

    void init(GridFunc<T>& gf_rhod)
    {
        PB<T>::epsilon_.Gepsilon_rho(gf_rhod);
        PB<T>::initialized_ = true;
    };

    void jacobi(GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W);

    void get_vepsilon(GridFunc<T>&, GridFunc<T>&, GridFunc<T>&);
};

} // namespace pb

#endif
