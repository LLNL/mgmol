// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef DELH4_H
#define DELH4_H

#include "FDoper.h"

namespace pb
{

template <class T>
class Delxh4 : public FDoper<T>
{
public:
    Delxh4(const Grid& mygrid) : FDoper<T>(mygrid) { }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B) override
    {
        this->del1_4th(A, B, 0);
    }

    ~Delxh4() override {};

    static short minNumberGhosts() { return 2; }
};
template <class T>
class Delyh4 : public FDoper<T>
{
public:
    Delyh4(const Grid& mygrid) : FDoper<T>(mygrid) { }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B) override
    {
        this->del1_4th(A, B, 1);
    }

    ~Delyh4() override {};

    static short minNumberGhosts() { return 2; }
};
template <class T>
class Delzh4 : public FDoper<T>
{
public:
    Delzh4(const Grid& mygrid) : FDoper<T>(mygrid) { }

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B) override
    {
        this->del1_4th(A, B, 2);
    }

    ~Delzh4() override {};

    static short minNumberGhosts() { return 2; }
};

} // namespace pb

#endif
