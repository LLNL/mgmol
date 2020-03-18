// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PBH2_H
#define PBH2_H

#include "PB.h"

// operator -epsilon div*grad -(grad epsilon)* grad O(h**2)

namespace pb
{
template <class T>
class PBh2 : public PB<T>
{

private:
    void pb_2nd(GridFunc<T>& A, GridFunc<T>& B);

public:
    // constructor
    PBh2(const Grid& mygrid, DielFunc<T>& myepsilon) : PB<T>(mygrid, myepsilon)
    {
        PB<T>::initialized_ = true;
    };
    PBh2(const Grid& mygrid, const double e0, const double rho0,
        const double drho0)
        : PB<T>(mygrid, e0, rho0, drho0)
    {
        PB<T>::initialized_ = false;
    };

    // construct a coarse grid operator
    PBh2 coarseOp(const Grid& mygrid);

    PBh2 replicatedOp(const Grid&);

    static short minNumberGhosts() { return 1; }

    void setLowerOrderGrid() override
    {
        FDoper<T>::setFDLowerOrderGrid(minNumberGhosts());
    }

    PBh2& getLowerOrderOp() { return *this; }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B) override
    {
        pb_2nd(A, B);
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
