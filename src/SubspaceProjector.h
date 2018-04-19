// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#ifndef SUBSPACEPROJ_H
#define SUBSPACEPROJ_H

#include "global.h"

class LocGridOrbitals;
class ProjectedMatricesInterface;
template <class T>
class SquareLocalMatrices;

class SubspaceProjector
{
private:
    LocGridOrbitals& subspace_;
    ProjectedMatricesInterface& proj_matrices_;
    short chromatic_number_;
    short subdivx_;
    int lda_;
    int loc_numpt_;
    
public:
    SubspaceProjector(LocGridOrbitals& subspace);

    ~SubspaceProjector(){}

    void projectOut(LocGridOrbitals&, SquareLocalMatrices<MATDTYPE>* mask=0);
};

#endif
