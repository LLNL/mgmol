// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$

class SymmetricPair
{
    int i1_;
    int i2_;

public:
    SymmetricPair()
    {
        i1_ = -1;
        i2_ = -1;
    }

    SymmetricPair(const int i1, const int i2) { setup(i1, i2); }

    void setup(const int i1, const int i2)
    {
        i1_ = (i1 <= i2) ? i1 : i2;
        i2_ = (i1_ == i1) ? i2 : i1;
    }

    bool operator<(const SymmetricPair& pair) const
    {
        if (i1_ < pair.i1_) return true;
        if (i1_ > pair.i1_) return false;
        if (i2_ < pair.i2_) return true;
        return false;
    }
};
