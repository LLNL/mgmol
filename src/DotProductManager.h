// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef MGMOL_DotProductManager_H
#define MGMOL_DotProductManager_H

template <class T>
class DotProductManager
{
public:
    DotProductManager() {};

    virtual ~DotProductManager() {};

    virtual double dotProduct(T& a, const T& b) = 0;
};

#endif
