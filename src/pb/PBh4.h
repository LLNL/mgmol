// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PB_H4_H
#define PB_H4_H

#include "PB.h"
#include "PBh2.h"

// operator -epsilon div*grad -(grad epsilon)* grad

namespace pb
{
template <class T>
class PBh4 : public PB<T>
{

private:
    void pb_4th(GridFunc<T>& A, GridFunc<T>& B);

    PBh2<T>* lower_order_op_;

public:
    // constructor
    PBh4(const Grid& mygrid, DielFunc<T>& myepsilon) : PB<T>(mygrid, myepsilon)
    {
        PB<T>::initialized_ = true;
        lower_order_op_     = nullptr;
    };
    PBh4(const Grid& mygrid, const double e0, const double rho0,
        const double drho0)
        : PB<T>(mygrid, e0, rho0, drho0)
    {
        PB<T>::initialized_ = false;
        lower_order_op_     = nullptr;
    };

    ~PBh4() override
    {
        if (lower_order_op_ != nullptr)
        {
            delete lower_order_op_;
            lower_order_op_ = nullptr;
        }
    }

    // construct a coarse grid operator
    PBh4 coarseOp(const Grid& mygrid)
    {
        // cout<<"construct a coarse grid operator"<<endl;
        Grid coarse_G = mygrid.coarse_grid();
        DielFunc<T> ecoarse(coarse_G, PB<T>::epsilon_.epsilon_max());
        PB<T>::epsilon_.restrict3D(ecoarse);

        PBh4 A(coarse_G, ecoarse);

        return A;
    }

    void setLowerOrderGrid() override
    {
        this->setFDLowerOrderGrid(PBh2<T>::minNumberGhosts());
    }

    PBh2<T>& getLowerOrderOp()
    {
        if (lower_order_op_ == nullptr)
        {
            this->setFDLowerOrderGrid(PBh2<T>::minNumberGhosts());
            lower_order_op_
                = new PBh2<T>(PB<T>::getLowerOrderGrid(), PB<T>::epsilon_);
        }
        return *lower_order_op_;
    }

    static short minNumberGhosts() { return 2; }

    PBh4 replicatedOp(const Grid&);

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B) override
    {
        pb_4th(A, B);
        B.set_bc(A.bc(0), A.bc(1), A.bc(2));
    }

    void init(GridFunc<T>& gf_rhod)
    {
        PB<T>::epsilon_.Gepsilon_rho(gf_rhod);
        PB<T>::initialized_ = true;
    };

    void jacobi(GridFunc<T>& A, const GridFunc<T>& B, GridFunc<T>& W) override;

    void get_vepsilon(GridFunc<T>&, GridFunc<T>&, GridFunc<T>&) override;
};

} // namespace pb

#endif
