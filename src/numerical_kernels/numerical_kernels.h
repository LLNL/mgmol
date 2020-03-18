// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <vector>

template <typename T1, typename T2, typename T3>
void nonOrthoRhoKernel(const short i0, const short ib, const short j0,
    const short jb, const int x0, const int xb, const T3* const mat,
    const int ld, const std::vector<const T1*>& psi, const double factor,
    T2* const rho);

template <typename T1, typename T2, typename T3>
void nonOrthoRhoKernelDiagonalBlock(const short i0, const short ib,
    const int x0, const int xb, const T3* const mat, const int ld,
    const std::vector<const T1*>& psi, T2* const rho);
