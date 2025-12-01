// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "OrthoAndersonMix.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "SquareLocalMatrices.h"

template <class T>
void OrthoAndersonMix<T>::postprocessUpdate()
{
    SquareLocalMatrices<double, MemorySpace::Host> ortho_transform(
        x_.subdivx(), x_.chromatic_number());
    x_.orthonormalizeLoewdin(false, &ortho_transform);
    for (int j = 0; j < m_; j++)
    {
        // xi_[j]->multiplyByMatrix(ortho_transform);
        // fi_[j]->multiplyByMatrix(ortho_transform);
    }
}

template class OrthoAndersonMix<ExtendedGridOrbitals<ORBDTYPE>>;
template class OrthoAndersonMix<LocGridOrbitals<ORBDTYPE>>;
