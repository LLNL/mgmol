// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef MGMOL_DotProductManagerWithDM_H
#define MGMOL_DotProductManagerWithDM_H

#include "DotProductManager.h"

template <class T>
class DotProductWithDM : public DotProductManager<T>
{
public:
    double dotProduct(T& phi0, const T& phi1) override;
};

#endif
