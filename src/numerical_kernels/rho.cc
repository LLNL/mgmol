// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "../global.h"
#include "Timer.h"
#include "numerical_kernels.h"

// #define WTIMERS

#ifdef WTIMERS
Timer nonOrthoRhoKernel_tm("nonOrthoRhoKernel");
Timer nonOrthoRhoKernelDiagonalBlock_tm("nonOrthoRhoKernelDiagonalBlock");
#endif

// Numerical kernel:
// Triple loop over pairs of functions and space
template <typename T1, typename T2, typename T3>
void nonOrthoRhoKernel(const short i0, const short ib, const short j0,
    const short jb, const int x0, const int xb, const T3* const mat,
    const int ld, const std::vector<const T1*>& psi, const double factor,
    T2* const rho)
{
#ifdef WTIMERS
    nonOrthoRhoKernel_tm.start();
#endif

    T2* __restrict__ prho = &rho[x0];

    // loop over function j
    for (short j = j0; j < j0 + jb; j++)
    {
        const T1* __restrict__ psij = &psi[j][x0];

        // loop over function i
        for (short i = i0; i < i0 + ib; i++)
        {
            const T1* __restrict__ psii = &psi[i][x0];

            const double kij = factor * (double)mat[j * ld + i];

            // loop over space
            for (int x = 0; x < xb; x++)
            {
                prho[x] += (T2)(kij * (double)psii[x] * (double)psij[x]);
            }
        }
    }

#ifdef WTIMERS
    nonOrthoRhoKernel_tm.stop();
#endif
}

template <typename T1, typename T2, typename T3>
void nonOrthoRhoKernelDiagonalBlock(const short i0, const short ib,
    const int x0, const int xb, const T3* const mat, const int ld,
    const std::vector<const T1*>& psi, T2* const rho)
{
#ifdef WTIMERS
    nonOrthoRhoKernelDiagonalBlock_tm.start();
#endif
    T2* __restrict__ prho = &rho[x0];

    // loop over function i
    for (short i = i0; i < i0 + ib; i++)
    {
        const T1* __restrict__ psii = &psi[i][x0];

        // loop over function j
        for (short j = i0; j < i; j++)
        {
            const T1* __restrict__ psij = &psi[j][x0];

            // factor 2 since we use symmetry and compute
            // only j<i elements
            const double kij = 2. * (double)mat[j * ld + i];

            // loop over space
            for (int x = 0; x < xb; x++)
            {
                prho[x] += (T2)(kij * (double)psii[x] * (double)psij[x]);
            }
        }

        const double kii = (double)mat[i * ld + i];
        // loop over space
        for (int x = 0; x < xb; x++)
        {
            prho[x] += (T2)(kii * (double)psii[x] * (double)psii[x]);
        }
    }

#ifdef WTIMERS
    nonOrthoRhoKernelDiagonalBlock_tm.stop();
#endif
}

template void nonOrthoRhoKernel(const short i0, const short ib, const short j0,
    const short jb, const int x0, const int xb, const MATDTYPE* const mat,
    const int ld, const std::vector<const ORBDTYPE*>& psi, const double factor,
    RHODTYPE* const rho);

template void nonOrthoRhoKernelDiagonalBlock(const short i0, const short ib,
    const int x0, const int xb, const MATDTYPE* const mat, const int ld,
    const std::vector<const ORBDTYPE*>& psi, RHODTYPE* const rho);
