// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef SQUARELOCALMATRICES_H
#define SQUARELOCALMATRICES_H

#include "LocalMatrices.h"

template <class T>
class SquareLocalMatrices : public LocalMatrices<T>
{
public:
    SquareLocalMatrices(const int subdiv, const int m);

    void fillUpperWithLower();
    void setDiagonal2Zero();
    void transpose();
};

#endif
