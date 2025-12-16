// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "DotProductWithInvS.h"

#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"
#include "ProjectedMatricesInterface.h"
#include "SquareLocalMatrices.h"

template <class T>
double DotProductWithInvS<T>::dotProduct(T& phi0, const T& phi1)
{
    Mesh* mymesh               = Mesh::instance();
    const int subdivx          = mymesh->subdivx();
    const int chromatic_number = phi0.chromatic_number();

    SquareLocalMatrices<MATDTYPE, MemorySpace::Host> ss(
        subdivx, chromatic_number);

    phi0.computeLocalProduct(phi1, ss);

    ProjectedMatricesInterface* proj_matrices = phi0.getProjMatrices();

    return proj_matrices->dotProductWithInvS(ss);
}

template class DotProductWithInvS<LocGridOrbitals<ORBDTYPE>>;
template class DotProductWithInvS<ExtendedGridOrbitals<ORBDTYPE>>;
