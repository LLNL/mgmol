// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#ifndef ORBITALS_H
#define ORBITALS_H

class Orbitals
{
    int iterative_index_;

public:
    Orbitals() { iterative_index_ = -10; }

    Orbitals(const Orbitals& A, const bool copy_data)
    {
        if (copy_data)
        {
            iterative_index_ = A.iterative_index_;
        }
        else
        {
            iterative_index_ = -10;
        }
    }

    void resetIterativeIndex() { iterative_index_ = 0; }

    void setIterativeIndex(const Orbitals& orbitals)
    {
        iterative_index_ = orbitals.iterative_index_;
    }

    void setIterativeIndex(const short iterative_index)
    {
#ifdef DEBUG
        MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
        mmpi.barrier();
#endif
        iterative_index_ = iterative_index;
    }

    void incrementIterativeIndex(const short inc = 1)
    {
        iterative_index_ += inc;
    };

    bool compareIterativeIndex(const Orbitals& orbitals) const
    {
        assert(orbitals.iterative_index_ >= 0);
        return (iterative_index_ == orbitals.iterative_index_);
    }

    short getIterativeIndex() const { return iterative_index_; }
};
#endif
