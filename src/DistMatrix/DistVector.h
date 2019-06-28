// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DISTVECTOR_H
#define MGMOL_DISTVECTOR_H

#include "DistMatrix.h"

#include <string>

namespace dist_matrix
{

template <class T>
class DistVector : public DistMatrix<T>
{
public:
    DistVector<T>(const std::string& name, const int m) :
        DistMatrix<T>(name, m, 1)
    {}

    void assign(const T* const v)
    {
        assignColumn(v, 0);
    }

};

}
#endif
