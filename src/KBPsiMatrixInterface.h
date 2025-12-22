// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_KBPSIMATRIX_INTERFACE_H
#define MGMOL_KBPSIMATRIX_INTERFACE_H

#include "Timer.h"
#include "VariableSizeMatrix.h"
#include "tools.h"

template <class T>
class ProjectedMatrices;
class ProjectedMatricesInterface;
class Ions;
class Ion;

class KBPsiMatrixInterface
{
    int iterative_index_;

protected:
    static Timer computeLocalElement_tm_;

public:
    KBPsiMatrixInterface() : iterative_index_(-1) {};

    virtual ~KBPsiMatrixInterface() {};

    int getIterativeIndex() const { return iterative_index_; }
    void setOutdated() { iterative_index_ = -1; }
    void setIterativeIndex(const int index)
    {
        assert(index >= 0);
#ifdef DEBUG
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        iterative_index_ = index;
    }

    virtual void addKBPsi(const int gid, const int st, const double val)  = 0;
    virtual void addKBBPsi(const int gid, const int st, const double val) = 0;

    virtual void printTimers(std::ostream& os);

    void computeLocalElement(Ion& ion, const int istate, const int iloc,
        const ORBDTYPE* const psi, const bool flag);
};

#endif
