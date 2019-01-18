// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_EIGENDMSTRATEGY_H
#define MGMOL_EIGENDMSTRATEGY_H

#include "DMStrategy.h"
class ProjectedMatricesInterface;

template <class T>
class EigenDMStrategy : public DMStrategy
{
private:
    T* current_orbitals_;
    ProjectedMatricesInterface* proj_matrices_;

public:
    EigenDMStrategy(T* current_orbitals,
        ProjectedMatricesInterface* proj_matrices);

    void initialize();
    int update();

    bool needH() const { return true; }

    void stripDM() {}
    void dressDM() {}
    void reset() {}
};

#endif
