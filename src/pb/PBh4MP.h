// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: PBh4MP.h,v 1.10 2010/01/28 22:56:47 jeanluc Exp $
#ifndef PB_MH4P_H
#define PB_MH4P_H

#include "PB.h"
#include "PBh4M.h"

namespace pb
{
template <class T>
class PBh4MP : public PBh4M<T>
{

public:
    // constructor
    PBh4MP(const Grid& mygrid, DielFunc<T>& myepsilon)
        : PBh4M<T>(mygrid, myepsilon) {};
    PBh4MP(const Grid& mygrid, const double e0, const double rho0,
        const double drho0)
        : PBh4M<T>(mygrid, e0, rho0, drho0) {};

    PBh4MP(const Grid& mygrid, DielFunc<T>& myepsilon, GridFunc<T>& pp)
        : PBh4M<T>(mygrid, myepsilon, pp) {};

    // construct a coarse grid operator
    PBh4MP coarseOp(const Grid& mygrid)
    {
        Grid coarse_G = mygrid.coarse_grid();
        DielFunc<T> ecoarse(coarse_G, PBh4M<T>::epsilon_.epsilon_max());
        GridFunc<T> ppcoarse(coarse_G, 1, 1, 1);
        PB<T>::epsilon_.restrict3D(ecoarse);
        PBh4M<T>::pp_.restrict3D(ppcoarse);

        PBh4MP A(coarse_G, ecoarse, ppcoarse);
        // A.epsilon_.set_grid(A.grid());

        return A;
    }

    PBh4MP replicatedOp(const Grid& replicated_grid)
    {
        T* replicated_func = new T[replicated_grid.gsize()];

        this->epsilon_.init_vect(replicated_func, 'g');

        DielFunc<T> replicated_epsilon(replicated_grid);
        replicated_epsilon.assign(replicated_func, 0);

        this->pp_.init_vect(replicated_func, 'g');
        GridFunc<T> replicated_pp(replicated_grid, 1, 1, 1);
        replicated_pp.assign(replicated_func);

        // cout<<"construct a coarse grid operator"<<endl;
        PBh4MP A(replicated_grid, replicated_epsilon, replicated_pp);
        // cout<<"New PBh4M coarse grid operator: My Grid has
        // "<<A.pp_.grid().sizeg()<<" points"<<endl;
        // A.epsilon_.set_grid(A.grid());

        return A;
    }

    void setLowerOrderGrid() override
    {
        FDoper<T>::setFDLowerOrderGrid(minNumberGhosts());
    }

    PBh4MP& lowerOrderOp() { return *this; }

    static short minNumberGhosts() { return 1; }

    void rhs(GridFunc<T>& A, GridFunc<T>& B) const override
    {
        FDoper<T>::rhs_4th_Mehr2(A, B);
        B.set_bc(A.bc(0), A.bc(1), A.bc(2));
    }
    void rhs(GridFunc<T>& A, T* const B) const override
    {
        FDoper<T>::rhs_4th_Mehr2(A, B);
    }
};

} // namespace pb

#endif
