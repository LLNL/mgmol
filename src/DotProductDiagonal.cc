// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "DotProductDiagonal.h"

#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"
#include "ProjectedMatricesInterface.h"
#include "SquareLocalMatrices.h"

template <>
double DotProductDiagonal<ExtendedGridOrbitals<ORBDTYPE>>::dotProduct(
    ExtendedGridOrbitals<ORBDTYPE>& phi0,
    const ExtendedGridOrbitals<ORBDTYPE>& phi1)
{
    const int chromatic_number = phi0.chromatic_number();
    std::vector<DISTMATDTYPE> ss(chromatic_number);
    phi0.computeDiagonalElementsDotProduct(phi1, ss);

    ProjectedMatricesInterface* proj_matrices = phi0.getProjMatrices();

    return proj_matrices->getTraceDiagProductWithInvS(ss);
}

template <>
double DotProductDiagonal<LocGridOrbitals<ORBDTYPE>>::dotProduct(
    LocGridOrbitals<ORBDTYPE>& phi0, const LocGridOrbitals<ORBDTYPE>& phi1)
{
    const int numst                           = phi0.numst();
    ProjectedMatricesInterface* proj_matrices = phi0.getProjMatrices();
    assert(proj_matrices != nullptr);

    std::vector<DISTMATDTYPE> ss;
    Control& ct = *(Control::instance());
    if (ct.short_sighted)
    {
        phi0.computeDiagonalElementsDotProductLocal(phi1, ss);
    }
    else
    {
        ss.resize(numst);
        phi0.computeDiagonalElementsDotProduct(phi1, ss);
    }

    return proj_matrices->getTraceDiagProductWithInvS(ss);
}
