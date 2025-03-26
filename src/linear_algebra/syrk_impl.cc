// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Timer.h"
#include "mputils.h"

#include "MGmol_blas1.h"
#include "syrk_impl.h"

#include "blas3_c.h"

#include <cassert>
#include <iostream>
#include <vector>

Timer mpsyrk_tm("mpsyrk");
Timer tttsyrk_tm("tttsyrk");

using LAU_H = LinearAlgebraUtils<MemorySpace::Host>;

template <>
void syrk_impl(const char uplo, const char trans, const int n, const int k,
    const double alpha, const float* const a, const int lda, const double beta,
    float* c, const int ldc)
{
    mpsyrk_tm.start();

    if (beta == 1. && (alpha == 0. || n == 0 || k == 0)) return;

    /* case Trans == 'N' */
    if (trans == 'N' || trans == 'n')
    {
        /* buffer to hold accumulation in double */
        std::vector<double> buff(n);
        if (uplo == 'U' || uplo == 'u')
        {
            for (int j = 0; j < n; j++)
            {
                const int len = j + 1;
                std::fill(buff.begin(), buff.begin() + len, 0.);
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const float* colL = a + lda * l;
                    /* get multiplier */
                    double mult = static_cast<double>(
                        alpha * colL[j]); // same as alpha * a[lda*l + j];
                    LAU_H::MPaxpy(len, mult, colL, buff.data());
                }
                /* Update col j of upper part of matrix C. */
                /* Get pointer to beginning of column j in C. */
                float* cj = c + ldc * j;
                LAU_H::MPscal(len, beta, cj);
                for (int i = 0; i < len; i++)
                    cj[i] += static_cast<float>(buff[i]);
            }
        }
        else /* uplo = 'L' or 'l' */
        {
            for (int j = 0; j < n; j++)
            {
                const int len = n - (j + 1);
                std::fill(buff.begin(), buff.begin() + len, 0.);
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const float* colL = a + lda * l + j;
                    /* get multiplier */
                    double mult = static_cast<double>(
                        alpha * colL[0]); // same as alpha * a[lda*l + j];
                    LAU_H::MPaxpy(len, mult, colL, buff.data());
                }
                /* Update col j of upper part of matrix C. */
                /* Get pointer to beginning of column j in C. */
                float* cj = c + ldc * j + j;
                LAU_H::MPscal(len, beta, cj);
                for (int i = 0; i < len; i++)
                    cj[i] += static_cast<float>(buff[i]);
            }
        }
    }
    else /* Trans == 'T' or 'C' */
    {
        if (uplo == 'U' || uplo == 'u')
        {
            for (int j = 0; j < n; j++)
            {
                const float* __restrict__ aj = a + lda * j;
                for (int i = 0; i < j; i++)
                {
                    const int pos                = ldc * j + i;
                    const float* __restrict__ ai = a + lda * i;
                    double bc = static_cast<double>(c[pos]) * beta;
                    c[pos]    = static_cast<float>(
                        alpha * LAU_H::MPdot(k, ai, aj) + bc);
                }
            }
        }
        else /* uplo = 'L' or 'l' */
        {
            for (int j = 0; j < n; j++)
            {
                const float* __restrict__ aj = a + lda * j;
                for (int i = j; i < n; i++)
                {
                    const int pos                = ldc * j + i;
                    const float* __restrict__ ai = a + lda * i;
                    double bc = static_cast<double>(c[pos]) * beta;
                    c[pos]    = static_cast<float>(
                        alpha * LAU_H::MPdot(k, ai, aj) + bc);
                }
            }
        }
    }
    mpsyrk_tm.stop();
}

template <typename T1, typename T2>
void syrk_impl(const char uplo, const char trans, const int n, const int k,
    const double alpha, const T1* const a, const int lda, const double beta,
    T2* c, const int ldc)
{
    tttsyrk_tm.start();

    if (beta == 1. && (alpha == 0. || n == 0 || k == 0)) return;

    /* case Trans == 'N' */
    if (trans == 'N' || trans == 'n')
    {
        /* buffer to hold accumulation in double */
        std::vector<double> buff(n);
        if (uplo == 'U' || uplo == 'u')
        {
            for (int j = 0; j < n; j++)
            {
                const int len = j + 1;
                std::fill(buff.begin(), buff.begin() + len, 0.);
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const T1* colL = a + lda * l;
                    /* get multiplier */
                    double mult = static_cast<double>(
                        alpha * colL[j]); // same as alpha * a[lda*l + j];
                    LAU_H::MPaxpy(len, mult, colL, buff.data());
                }
                /* Update col j of upper part of matrix C. */
                /* Get pointer to beginning of column j in C. */
                T2* cj = c + ldc * j;
                LAU_H::MPscal(len, beta, cj);
                for (int i = 0; i < len; i++)
                    cj[i] += (T2)buff[i];
            }
        }
        else /* uplo = 'L' or 'l' */
        {
            for (int j = 0; j < n; j++)
            {
                const int len = n - (j + 1);
                std::fill(buff.begin(), buff.begin() + len, 0.);
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const T1* colL = a + lda * l + j;
                    /* get multiplier */
                    double mult = static_cast<double>(
                        alpha * colL[0]); // same as alpha * a[lda*l + j];
                    LAU_H::MPaxpy(len, mult, colL, buff.data());
                }
                /* Update col j of upper part of matrix C. */
                /* Get pointer to beginning of column j in C. */
                T2* cj = c + ldc * j + j;
                LAU_H::MPscal(len, beta, cj);
                for (int i = 0; i < len; i++)
                    cj[i] += (T2)buff[i];
            }
        }
    }
    else /* Trans == 'T' or 'C' */
    {
        if (uplo == 'U' || uplo == 'u')
        {
            for (int j = 0; j < n; j++)
            {
                const T1* __restrict__ aj = a + lda * j;
                for (int i = 0; i < j; i++)
                {
                    const int pos             = ldc * j + i;
                    const T1* __restrict__ ai = a + lda * i;
                    double bc = static_cast<double>(c[pos]) * beta;
                    c[pos]
                        = static_cast<T2>(alpha * LAU_H::MPdot(k, ai, aj) + bc);
                }
            }
        }
        else /* uplo = 'L' or 'l' */
        {
            for (int j = 0; j < n; j++)
            {
                const T1* __restrict__ aj = a + lda * j;
                for (int i = j; i < n; i++)
                {
                    const int pos             = ldc * j + i;
                    const T1* __restrict__ ai = a + lda * i;
                    double bc = static_cast<double>(c[pos]) * beta;
                    c[pos]
                        = static_cast<T2>(alpha * LAU_H::MPdot(k, ai, aj) + bc);
                }
            }
        }
    }

    tttsyrk_tm.stop();
}

template void syrk_impl<double, float>(const char uplo, const char trans,
    const int n, const int k, const double alpha, const double* const a,
    const int lda, const double beta, float* c, const int ldc);
template void syrk_impl<float, double>(const char uplo, const char trans,
    const int n, const int k, const double alpha, const float* const a,
    const int lda, const double beta, double* c, const int ldc);
