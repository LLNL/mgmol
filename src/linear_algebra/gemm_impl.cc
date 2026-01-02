// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "gemm_impl.h"
#include "Timer.h"
#include "mputils.h"

#include <cassert>
#include <vector>

Timer mpgemm_tm("mpgemm");
Timer tttgemm_tm("tttgemm");

using LAU_H = LinearAlgebraUtils<MemorySpace::Host>;

template <typename T1, typename T2, typename T3>
void gemm_impl(const char transa, const char transb, const int m, const int n,
    const int k, const double alpha, const T1* const a, const int lda,
    const T2* const b, const int ldb, const double beta, T3* const c,
    const int ldc)
{
    tttgemm_tm.start();
    // std::cout<<"template MPgemm..."<<std::endl;

    if (beta == 1. && (alpha == 0. || m == 0 || n == 0 || k == 0)) return;

    /* case transb == 'N' and transa == 'N' */
    if (transb == 'N' || transb == 'n')
    {
        if (transa == 'N' || transa == 'n')
        {
            /* buffer to hold accumulation in double */
            std::vector<double> buff(m);
            for (int j = 0; j < n; j++)
            {
                std::fill(buff.begin(), buff.end(), 0.);
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const T1* colL = a + lda * l;
                    /* get multiplier */
                    double mult = (double)(alpha * b[ldb * j + l]);
                    LAU_H::MPaxpy(m, mult, colL, buff.data());
                }
                /* Update col j of of result matrix C. */
                /* Get pointer to beginning of column j in C. */
                T3* cj = c + ldc * j;
                LAU_H::MPscal(m, beta, cj);
                for (int i = 0; i < m; i++)
                {
                    cj[i] += (T3)buff[i];
                }
            }
        }
        else /* transa == 'T'/'C' */
        {
            for (int j = 0; j < n; j++)
            {
                const T2* __restrict__ bj = b + ldb * j;
                for (int i = 0; i < m; i++)
                {
                    const int pos = ldc * j + i;
                    double bc     = static_cast<double>(c[pos]) * beta;
                    const T1* __restrict__ ai = a + lda * i;
                    c[pos]
                        = static_cast<T3>(alpha * LAU_H::MPdot(k, ai, bj) + bc);
                }
            }
        }
    }
    else /* transb == 'T'/'C' */
    {
        if (transa == 'N' || transa == 'n')
        {
            /* buffer to hold accumulation in double */
            std::vector<double> buff(m);
            for (int j = 0; j < n; j++)
            {
                std::fill(buff.begin(), buff.end(), 0.);
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const T1* colL = a + lda * l;
                    /* get multiplier */
                    double mult = (double)(alpha * b[ldb * l + j]);
                    LAU_H::MPaxpy(m, mult, colL, buff.data());
                }
                /* Update col j of of result matrix C. */
                /* Get pointer to beginning of column j in C. */
                T3* cj = c + ldc * j;
                LAU_H::MPscal(m, beta, cj);
                for (int i = 0; i < m; i++)
                {
                    cj[i] += (T3)buff[i];
                }
            }
        }
        else /* transa == 'T'/'C' */
        {
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    const int pos = ldc * j + i;
                    const T1* ai  = a + lda * i;
                    double sum    = 0.;
                    for (int l = 0; l < k; l++)
                    {
                        sum += alpha * ai[l] * b[ldb * l + j];
                    }
                    sum += (double)(beta * c[pos]);
                    c[pos] = (T3)sum;
                }
            }
        }
    }

    tttgemm_tm.stop();
}

// input/output in float, computation in double
template <>
void gemm_impl<float, float, float>(const char transa, const char transb,
    const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const float* const b, const int ldb,
    const double beta, float* const c, const int ldc)
{
    mpgemm_tm.start();

    if (beta == 1. && (alpha == 0. || m == 0 || n == 0 || k == 0)) return;

    /* case transb == 'N' and transa == 'N' */
    if (transb == 'N' || transb == 'n')
    {
        if (transa == 'N' || transa == 'n')
        {
            /* buffer to hold accumulation in double */
            std::vector<double> buff(m);
            for (int j = 0; j < n; j++)
            {
                std::fill(buff.begin(), buff.end(), 0);
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const float* colL = a + lda * l;
                    /* get multiplier */
                    double mult = (double)(alpha * b[ldb * j + l]);
                    LAU_H::MPaxpy(m, mult, colL, buff.data());
                }
                /* Update col j of of result matrix C. */
                /* Get pointer to beginning of column j in C. */
                float* cj = c + ldc * j;
                LAU_H::MPscal(m, beta, cj);
                for (int i = 0; i < m; i++)
                    cj[i] += (float)buff[i];
            }
        }
        else /* transa == 'T'/'C' */
        {
            for (int j = 0; j < n; j++)
            {
                const float* __restrict__ bj = b + ldb * j;
                for (int i = 0; i < m; i++)
                {
                    const int pos                = ldc * j + i;
                    double bc                    = (double)c[pos] * beta;
                    const float* __restrict__ ai = a + lda * i;
                    c[pos] = (float)(alpha * LAU_H::MPdot(k, ai, bj) + bc);
                }
            }
        }
    }
    else /* transb == 'T'/'C' */
    {
        if (transa == 'N' || transa == 'n')
        {
            /* buffer to hold accumulation in double */
            std::vector<double> buff(m);
            for (int j = 0; j < n; j++)
            {
                std::fill(buff.begin(), buff.end(), 0);
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const float* colL = a + lda * l;
                    /* get multiplier */
                    double mult = (double)(alpha * b[ldb * l + j]);
                    LAU_H::MPaxpy(m, mult, colL, buff.data());
                }
                /* Update col j of of result matrix C. */
                /* Get pointer to beginning of column j in C. */
                float* cj = c + ldc * j;
                LAU_H::MPscal(m, beta, cj);
                for (int i = 0; i < m; i++)
                    cj[i] += (float)buff[i];
            }
        }
        else /* transa == 'T'/'C' */
        {
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    const int pos   = ldc * j + i;
                    const float* ai = a + lda * i;
                    double sum      = 0.;
                    for (int l = 0; l < k; l++)
                    {
                        sum += alpha * ai[l] * b[ldb * l + j];
                    }
                    sum += (double)(beta * c[pos]);
                    c[pos] = (float)sum;
                }
            }
        }
    }

    mpgemm_tm.stop();
}

template void gemm_impl<double, float, double>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const double* const a, const int lda,
    const float* const b, const int ldb, const double beta, double* const c,
    const int ldc);
template void gemm_impl<float, double, float>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const float* const a, const int lda,
    const double* const b, const int ldb, const double beta, float* const c,
    const int ldc);
template void gemm_impl<double, double, float>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const double* const a, const int lda,
    const double* const b, const int ldb, const double beta, float* const c,
    const int ldc);
template void gemm_impl<float, float, double>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const float* const a, const int lda,
    const float* const b, const int ldb, const double beta, double* const c,
    const int ldc);
