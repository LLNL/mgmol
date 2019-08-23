// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_NONORTHODMSTRATEGY_H
#define MGMOL_NONORTHODMSTRATEGY_H

#include "DMStrategy.h"
#include "ProjectedMatricesInterface.h"

template <class T>
class NonOrthoDMStrategy : public DMStrategy
{
private:
    T* orbitals_;
    ProjectedMatricesInterface* proj_matrices_;
    const double mix_;

public:
    NonOrthoDMStrategy(T* orbitals, ProjectedMatricesInterface* proj_matrices,
        const double mix);

    void initialize();
    int update();

    bool needH() const { return true; }

    void stripDM();

    void dressDM();

    void reset() { proj_matrices_->resetDM(); }
};

#endif
