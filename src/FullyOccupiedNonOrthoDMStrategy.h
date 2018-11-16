// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef FULLYOCCUPIEDNONORTHODMSTRATEGY_H
#define FULLYOCCUPIEDNONORTHODMSTRATEGY_H

#include "DMStrategy.h"
class ProjectedMatricesInterface;

class FullyOccupiedNonOrthoDMStrategy : public DMStrategy
{
    ProjectedMatricesInterface* proj_matrices_;

public:
    FullyOccupiedNonOrthoDMStrategy(ProjectedMatricesInterface* proj_matrices);

    void initialize();
    int update();

    bool needH() const { return false; }

    void stripDM() {}
    void dressDM() {}
    void reset() {}
};

#endif
