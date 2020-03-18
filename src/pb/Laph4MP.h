// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PB_LAPH4MP_H
#define PB_LAPH4MP_H

#include "Laph4M.h"
namespace pb
{
template <class T>
class Laph4MP : public Laph4M<T>
{

public:
    Laph4MP(const Grid& mygrid) : Laph4M<T>(mygrid)
    {
        Lap<T>::name_ = "SPD Mehrstellen 4th order";
    };

    // construct a coarse grid operator
    Laph4MP coarseOp(const Grid& mygrid)
    {
        Grid coarse_G(mygrid.coarse_grid(), 1);

        Laph4MP A(coarse_G);

        return A;
    }

    Laph4MP replicatedOp(const Grid& replicated_grid)
    {
        Laph4MP replicated_A(replicated_grid);

        return replicated_A;
    }

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
