// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "mputils.h"
#include "MPIdata.h"
#include "Timer.h"
#include "magma_singleton.h"

#ifdef HAVE_MAGMA
#include <magma_v2.h>
#endif

#include <cassert>
#include <iostream>
#include <vector>

Timer dgemm_tm("dgemm");
Timer sgemm_tm("sgemm");
Timer mpgemm_tm("mpgemm");
Timer tttgemm_tm("tttgemm");

Timer dsyrk_tm("dsyrk");
Timer ssyrk_tm("ssyrk");
Timer mpsyrk_tm("mpsyrk");
Timer tttsyrk_tm("tttsyrk");

Timer mpdot_tm("mpdot");
Timer ttdot_tm("ttdot");

/* Function definitions. See mputils.h for comments */

using LAU_H = LinearAlgebraUtils<MemorySpace::Host>;
template <typename ScalarType>
using MemoryH = MemorySpace::Memory<ScalarType, MemorySpace::Host>;

#ifdef HAVE_MAGMA
using LAU_D = LinearAlgebraUtils<MemorySpace::Device>;
template <typename ScalarType>
using MemoryD = MemorySpace::Memory<ScalarType, MemorySpace::Device>;
#endif

/////////////////////////////
//          MPscal         //
/////////////////////////////
// MemorySpace::Host
template <>
void LAU_H::MPscal(const int len, const double scal, double* dptr)
{
#ifdef HAVE_MAGMA
    assert(magma_is_devptr(dptr) == 0);
#endif

    const int one = 1;
    DSCAL(&len, &scal, dptr, &one);
}

template <>
void LAU_H::MPscal(const int len, const double scal, float* dptr)
{
#ifdef HAVE_MAGMA
    assert(magma_is_devptr(dptr) == 0);
#endif

    if (scal == 1.)
        return;
    else if (scal == 0.)
    {
        memset(dptr, 0, len * sizeof(float));
    }
    else
    {
        for (int k = 0; k < len; k++)
        {
            double val = static_cast<double>(dptr[k]);
            dptr[k]    = static_cast<float>(scal * val);
        }
    }
}

// MemorySpace::Device
#ifdef HAVE_MAGMA
template <>
void LAU_D::MPscal(const int len, const double scal, double* dptr)
{
    assert(magma_is_devptr(dptr) == 1);

    int const increment   = 1;
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();
    magma_dscal(len, scal, dptr, increment, magma_singleton.queue_);
}

template <>
void LAU_D::MPscal(const int len, const double scal, float* dptr)
{
    assert(magma_is_devptr(dptr) == 1);

    if (scal == 1.)
        return;
    else if (scal == 0)
    {
        MemoryD<float>::set(dptr, 0, len);
    }
    else
    {
#ifdef HAVE_OPENMP_OFFLOAD
        float* dptr_alias = dptr;
#else
        std::unique_ptr<float[], void (*)(float*)> dptr_alias(
            MemoryH<float>::allocate(len), MemoryH<float>::free);
        MemorySpace::copy_to_host(dptr, dptr_alias);
#endif

        MGMOL_PARALLEL_FOR(dptr_alias)
        for (int k = 0; k < len; k++)
        {
            double val    = static_cast<double>(dptr_alias[k]);
            dptr_alias[k] = static_cast<float>(scal * val);
        }

#ifndef HAVE_OPENMP_OFFLOAD
        MemorySpace::copy_to_dev(dptr_alias, len, dptr);
#endif
    }
}
#endif

////////////////////////////
//          MPdot         //
////////////////////////////
// MemorySpace::Host
template <>
template <>
double LAU_H::MPdot<double, double>(
    const int len, const double* const xptr, const double* const yptr)
{
#ifdef HAVE_MAGMA
    assert(magma_is_devptr(xptr) == 0);
    assert(magma_is_devptr(yptr) == 0);
#endif

    const int one = 1;
    return DDOT(&len, xptr, &one, yptr, &one);
}

template <>
template <typename T1, typename T2>
double LAU_H::MPdot(
    const int len, const T1* __restrict__ xptr, const T2* __restrict__ yptr)
{
#ifdef HAVE_MAGMA
    assert(magma_is_devptr(xptr) == 0);
    assert(magma_is_devptr(yptr) == 0);
#endif

    mpdot_tm.start();

    double dot = 0.;
    for (int k = 0; k < len; k++)
    {
        double val1 = static_cast<double>(xptr[k]);
        double val2 = static_cast<double>(yptr[k]);
        dot += val1 * val2;
    }

    mpdot_tm.stop();

    return dot;
}

// MemorySpace::Device
#ifdef HAVE_MAGMA
template <>
template <>
double LAU_D::MPdot<double, double>(
    const int len, const double* const xptr, const double* const yptr)
{
    assert(magma_is_devptr(xptr) == 1);
    assert(magma_is_devptr(yptr) == 1);

    const int increment   = 1;
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();
    return magma_ddot(
        len, xptr, increment, yptr, increment, magma_singleton.queue_);
}

template <>
template <typename T1, typename T2>
double LAU_D::MPdot(
    const int len, const T1* __restrict__ xptr, const T2* __restrict__ yptr)
{
    assert(magma_is_devptr(xptr) == 1);
    assert(magma_is_devptr(yptr) == 1);

#ifndef HAVE_OPENMP_OFFLOAD
    std::unique_ptr<T1[], void (*)(T1*)> xptr_host(
        MemoryH<T1>::allocate(len), MemoryH<T1>::free);
    std::unique_ptr<T2[], void (*)(T2*)> yptr_host(
        MemoryH<T2>::allocate(len), MemoryH<T2>::free);
    copy_to_host(xptr, xptr_host, len);
    copy_to_host(yptr, yptr_host, len);
    return LAU_H::MPdot(len, x_ptr_host.get(), y_ptr_host.get());
#else
    double dot = 0.;
    // clang-format off
#pragma omp target teams distribute parallel for map(tofrom: dot) is_device_ptr(xptr, yptr)
    // clang-format on
    for (int k = 0; k < len; k++)
    {
        double val1 = static_cast<double>(xptr[k]);
        double val2 = static_cast<double>(yptr[k]);
        dot += val1 * val2;
    }

    return dot;
#endif
}
#endif

///////////////////////////////
////          MPaxpy         //
///////////////////////////////
// MemorySpace::Host
template <>
void LAU_H::MPaxpy(const int len, double scal, const double* __restrict__ xptr,
    double* __restrict__ yptr)
{
#ifdef HAVE_MAGMA
    assert(magma_is_devptr(xptr) == 0);
    assert(magma_is_devptr(yptr) == 0);
#endif

    const int one = 1;
    DAXPY(&len, &scal, xptr, &one, yptr, &one);
}

template <>
template <typename T1, typename T2>
void LAU_H::MPaxpy(const int len, double scal, const T1* __restrict__ xptr,
    T2* __restrict__ yptr)
{
#ifdef HAVE_MAGMA
    assert(magma_is_devptr(xptr) == 0);
    assert(magma_is_devptr(yptr) == 0);
#endif

#pragma omp parallel for simd
    for (int k = 0; k < len; k++)
    {
        yptr[k] += static_cast<T2>(scal * static_cast<double>(xptr[k]));
    }
}

// MemorySpace::Device
#ifdef HAVE_MAGMA
template <>
void LAU_D::MPaxpy(const int len, double scal, const double* __restrict__ xptr,
    double* __restrict__ yptr)
{
    assert(magma_is_devptr(xptr) == 1);
    assert(magma_is_devptr(yptr) == 1);

    const int increment   = 1;
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();
    return magma_daxpy(
        len, scal, xptr, increment, yptr, increment, magma_singleton.queue_);
}

template <>
template <typename T1, typename T2>
void LAU_D::MPaxpy(const int len, double scal, const T1* __restrict__ xptr,
    T2* __restrict__ yptr)
{
    assert(magma_is_devptr(xptr) == 1);
    assert(magma_is_devptr(yptr) == 1);

#ifndef HAVE_OPENMP_OFFLOAD
    std::unique_ptr<T1[], void (*)(double*)> xptr_host(
        MemoryH<T1>::allocate(len), MemoryH<T1>::free);
    std::unique_ptr<T2[], void (*)(float*)> yptr_host(
        MemoryH<T2>::allocate(len), MemoryH<T2>::free);
    copy_to_host(xptr, xptr_host, len);
    copy_to_host(yptr, yptr_host, len);
    LAU_H::MPaxpy(len, scal, xptr_host.get(), yptr_host.get());
    copy_to_dev(yptr_host, yptr, len);
#else
    // clang-format off
#pragma omp target teams distribute parallel for map(to: scal) is_device_ptr(xptr, yptr)
    // clang-format on
    for (int k = 0; k < len; k++)
    {
        yptr[k] += static_cast<T2>(scal * static_cast<double>(xptr[k]));
    }
#endif
}
#endif

void Tscal(const int len, const double scal, double* dptr)
{
    const int one = 1;
    DSCAL(&len, &scal, dptr, &one);
}
void Tscal(const int len, const float scal, float* dptr)
{
    const int one = 1;
    SSCAL(&len, &scal, dptr, &one);
}

double Tnrm2(const int len, const double* const dptr)
{
    const int one = 1;
    double nrm;
    nrm = DNRM2(&len, dptr, &one);

    return nrm;
}
float Tnrm2(const int len, const float* const dptr)
{
    const int one = 1;
    float nrm;
    nrm = SNRM2(&len, dptr, &one);

    return nrm;
}
double MPnrm2(const int len, const float* const dptr)
{
    double nrm = 0.;
    for (int k = 0; k < len; k++)
    {
        double val = (double)dptr[k];
        nrm += val * val;
    }
    return sqrt(nrm);
}

double Tdot(const int len, const double* const xptr, const double* const yptr)
{
    const int one = 1;
    double dot;
    dot = DDOT(&len, xptr, &one, yptr, &one);

    return dot;
}
float Tdot(const int len, const float* const xptr, const float* const yptr)
{
    const int one = 1;
    float dot;
    dot = SDOT(&len, xptr, &one, yptr, &one);

    return dot;
}

void Taxpy(const int len, double scal, const double* const xptr, double* yptr)
{
    const int one = 1;
    DAXPY(&len, &scal, xptr, &one, yptr, &one);
}
void Taxpy(const int len, float scal, const float* const xptr, float* yptr)
{
    const int one = 1;
    SAXPY(&len, &scal, xptr, &one, yptr, &one);
}

void Ttrsm(const char side, const char uplo, const char transa, const char diag,
    const int m, const int n, const double alpha, const double* const a,
    const int lda, double* const b, const int ldb)
{
    DTRSM(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

void Ttrsm(const char side, const char uplo, const char transa, const char diag,
    const int m, const int n, const float alpha, const float* const a,
    const int lda, float* const b, const int ldb)
{
    STRSM(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

void Tsyrk(const char uplo, const char trans, const int n, const int k,
    const double alpha, const double* const a, const int lda, const double beta,
    double* c, const int ldc)
{
    dsyrk_tm.start();
    DSYRK(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
    dsyrk_tm.stop();
}
void Tsyrk(const char uplo, const char trans, const int n, const int k,
    const float alpha, const float* const a, const int lda, const float beta,
    float* c, const int ldc)
{
    ssyrk_tm.start();
    SSYRK(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
    ssyrk_tm.stop();
}

void MPsyrk(const char uplo, const char trans, const int n, const int k,
    const double alpha, const double* const a, const int lda, const double beta,
    double* c, const int ldc)
{
    dsyrk_tm.start();
    DSYRK(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
    dsyrk_tm.stop();
}

void MPsyrk(const char uplo, const char trans, const int n, const int k,
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
                    double mult
                        = (double)(alpha
                                   * colL[j]); // same as alpha * a[lda*l + j];
                    LAU_H::MPaxpy(len, mult, colL, buff.data());
                }
                /* Update col j of upper part of matrix C. */
                /* Get pointer to beginning of column j in C. */
                float* cj = c + ldc * j;
                LAU_H::MPscal(len, beta, cj);
                for (int i = 0; i < len; i++)
                    cj[i] += (float)buff[i];
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
                    double mult
                        = (double)(alpha
                                   * colL[0]); // same as alpha * a[lda*l + j];
                    LAU_H::MPaxpy(len, mult, colL, buff.data());
                }
                /* Update col j of upper part of matrix C. */
                /* Get pointer to beginning of column j in C. */
                float* cj = c + ldc * j + j;
                LAU_H::MPscal(len, beta, cj);
                for (int i = 0; i < len; i++)
                    cj[i] += (float)buff[i];
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

void Tgemv(const char trans, const int m, const int n, const double alpha,
    const double* const a, const int lda, const double* const x, const int incx,
    const double beta, double* const y, const int incy)
{
    DGEMV(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}
void Tgemv(const char trans, const int m, const int n, const float alpha,
    const float* const a, const int lda, const float* const x, const int incx,
    const float beta, float* const y, const int incy)
{
    SGEMV(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

/////////////////////////////
//          MPgemm         //
/////////////////////////////

// MemorySpaceType
template <typename MemorySpaceType>
template <typename T1, typename T2, typename T3>
void LinearAlgebraUtils<MemorySpaceType>::MPgemm(const char /*transa*/,
    const char /*transb*/, const int /*m*/, const int /*n*/, const int /*k*/,
    const double /*alpha*/, const T1* const /*a*/, const int /*lda*/,
    const T2* const /*b*/, const int /*ldb*/, const double /*beta*/,
    T3* const /*c*/, const int /*ldc*/)
{
    assert(false);
}

// MemorySpace::Host
template <>
template <typename T1, typename T2, typename T3>
void LAU_H::MPgemm(const char transa, const char transb, const int m,
    const int n, const int k, const double alpha, const T1* const a,
    const int lda, const T2* const b, const int ldb, const double beta,
    T3* const c, const int ldc)
{
#ifdef HAVE_MAGMA
    assert(magma_is_devptr(a) == 0);
    assert(magma_is_devptr(b) == 0);
    assert(magma_is_devptr(c) == 0);
#endif

    tttgemm_tm.start();
    // if(onpe0)cout<<"template MPgemm..."<<endl;

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

// input/output in double, computation in double
template <>
template <>
void LAU_H::MPgemm<double, double, double>(const char transa, const char transb,
    const int m, const int n, const int k, const double alpha,
    const double* const a, const int lda, const double* const b, const int ldb,
    const double beta, double* const c, const int ldc)
{
#ifdef HAVE_MAGMA
    assert(magma_is_devptr(a) == 0);
    assert(magma_is_devptr(b) == 0);
    assert(magma_is_devptr(c) == 0);
#endif

    dgemm_tm.start();
    DGEMM(
        &transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    dgemm_tm.stop();
}

// input/output in float, computation in double
template <>
template <>
void LAU_H::MPgemm<float, float, float>(const char transa, const char transb,
    const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const float* const b, const int ldb,
    const double beta, float* const c, const int ldc)
{
#ifdef HAVE_MAGMA
    assert(magma_is_devptr(a) == 0);
    assert(magma_is_devptr(b) == 0);
    assert(magma_is_devptr(c) == 0);
#endif

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
                    c[pos] = (float)(alpha * MPdot(k, ai, bj) + bc);
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

// MemorySpace::Device
#ifdef HAVE_MAGMA
template <>
template <typename T1, typename T2, typename T3>
void LAU_D::MPgemm(const char transa, const char transb, const int m,
    const int n, const int k, const double alpha, const T1* const a,
    const int lda, const T2* const b, const int ldb, const double beta,
    T3* const c, const int ldc)
{
    assert(magma_is_devptr(a) == 1);
    assert(magma_is_devptr(b) == 1);
    assert(magma_is_devptr(c) == 1);

    std::vector<T1> a_host(lda * k);
    std::vector<T2> b_host(ldb * n);
    std::vector<T3> c_host(ldc * n);

    // Move the data to the host
    MemorySpace::copy_to_host(a, a_host);
    MemorySpace::copy_to_host(b, b_host);

    LAU_H::MPgemm(transa, transb, m, n, k, alpha, a_host.data(), lda,
        b_host.data(), ldb, beta, c_host.data(), ldc);

    // Move the data to the device
    MemorySpace::copy_to_dev(c_host, c);
}

// input/output in double, computation in double
template <>
template <>
void LAU_D::MPgemm(const char transa, const char transb, const int m,
    const int n, const int k, const double alpha, const double* const a,
    const int lda, const double* const b, const int ldb, const double beta,
    double* const c, const int ldc)
{
    assert(magma_is_devptr(a) == 1);
    assert(magma_is_devptr(b) == 1);
    assert(magma_is_devptr(c) == 1);

    dgemm_tm.start();
    // Transform char to magma_trans_t
    auto convert_to_magma_trans = [](const char trans) {
        if ((trans == 'N') || trans == 'n')
            return MagmaNoTrans;
        else if ((trans == 'T') || trans == 't')
            return MagmaTrans;
        else if ((trans == 'C') || trans == 'c')
            return MagmaConjTrans;
        else
        {
            std::cerr << "Unknown tranpose operation: " << trans << std::endl;
            return MagmaNoTrans;
        }
    };

    magma_trans_t magma_transa = convert_to_magma_trans(transa);
    magma_trans_t magma_transb = convert_to_magma_trans(transb);

    // Perform dgemm
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();
    magmablas_dgemm(magma_transa, magma_transb, m, n, k, alpha, a, lda, b, ldb,
        beta, c, ldc, magma_singleton.queue_);
    dgemm_tm.stop();
}
#endif

///////////////////////////////
//          MPgemmNN         //
///////////////////////////////

template <typename MemorySpaceType>
template <typename T1, typename T2, typename T3>
void LinearAlgebraUtils<MemorySpaceType>::MPgemmNN(const int m, const int n,
    const int k, const double alpha, const T1* const a, const int lda,
    const T2* const b, const int ldb, const double beta, T3* const c,
    const int ldc)
{
    char transa = 'n';
    char transb = 'n';
    LinearAlgebraUtils<MemorySpaceType>::MPgemm(
        transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

// input in float, computation in double
void MPgemmTN(const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const float* const b, const int ldb,
    const double beta, double* const c, const int ldc)
{
    char transa = 't';
    char transb = 'n';

    LAU_H::MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

// input in float, computation in double
void MPgemmTN(const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const float* const b, const int ldb,
    const double beta, float* const c, const int ldc)
{
    char transa = 't';
    char transb = 'n';

    LAU_H::MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

/////// additional calls ... may be removed later if unused

template <typename T1, typename T2, typename T3>
void MPgemmTN(const int m, const int n, const int k, const double alpha,
    const T1* const a, const int lda, const T2* const b, const int ldb,
    const double beta, T3* const c, const int ldc)
{
    // if(onpe0)cout<<"template MPgemmNN..."<<endl;
    char transa = 't';
    char transb = 'n';
    LAU_H::MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

template <typename T1, typename T2>
void MPsyrk(const char uplo, const char trans, const int n, const int k,
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
                    double mult
                        = (double)(alpha
                                   * colL[j]); // same as alpha * a[lda*l + j];
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
                    double mult
                        = (double)(alpha
                                   * colL[0]); // same as alpha * a[lda*l + j];
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

void MPcpy(double* const dest, const double* const src, const int n)
{
    memcpy(dest, src, n * sizeof(double));
}
void MPcpy(float* const dest, const float* const src, const int n)
{
    memcpy(dest, src, n * sizeof(float));
}
void MPcpy(
    double* __restrict__ dest, const float* __restrict__ src, const int n)
{
    for (int i = 0; i < n; i++)
        dest[i] = src[i];
}
void MPcpy(
    float* __restrict__ dest, const double* __restrict__ src, const int n)
{
    for (int i = 0; i < n; i++)
        dest[i] = src[i];
}

template void MPsyrk<double, float>(const char uplo, const char trans,
    const int n, const int k, const double alpha, const double* const a,
    const int lda, const double beta, float* c, const int ldc);
template void MPsyrk<float, double>(const char uplo, const char trans,
    const int n, const int k, const double alpha, const float* const a,
    const int lda, const double beta, double* c, const int ldc);

template void LAU_H::MPgemm<double, float, double>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const double* const a, const int lda,
    const float* const b, const int ldb, const double beta, double* const c,
    const int ldc);
template void LAU_H::MPgemm<float, double, float>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const float* const a, const int lda,
    const double* const b, const int ldb, const double beta, float* const c,
    const int ldc);
template void LAU_H::MPgemm<double, double, float>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const double* const a, const int lda,
    const double* const b, const int ldb, const double beta, float* const c,
    const int ldc);
template void LAU_H::MPgemm<float, float, double>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const float* const a, const int lda,
    const float* const b, const int ldb, const double beta, double* const c,
    const int ldc);
template void LAU_H::MPgemmNN<float, double, float>(const int m, const int n,
    const int k, const double alpha, const float* const a, const int lda,
    const double* const b, const int ldb, const double beta, float* const c,
    const int ldc);
template void LAU_H::MPgemmNN<double, double, double>(const int m, const int n,
    const int k, const double alpha, const double* const a, const int lda,
    const double* const b, const int ldb, const double beta, double* const c,
    const int ldc);
template void LAU_H::MPscal(const int len, const double scal, double* ptr);
template void LAU_H::MPscal(const int len, const double scal, float* ptr);
template double LAU_H::MPdot<double, double>(
    const int len, const double* const xptr, const double* const yptr);
template double LAU_H::MPdot<float, float>(const int len,
    const float* __restrict__ xptr, const float* __restrict__ yptr);
template double LAU_H::MPdot<double, float>(
    const int len, const double* const xptr, const float* const yptr);
template double LAU_H::MPdot<float, double>(
    const int len, const float* const xptr, const double* const yptr);
template void LAU_H::MPaxpy<float, double>(const int len, const double scal,
    const float* __restrict__ xptr, double* __restrict__ yptr);
template void LAU_H::MPaxpy<float, float>(const int len, const double scal,
    const float* __restrict__ xptr, float* __restrict__ yptr);

template void MPgemmTN<double, double, double>(const int m, const int n,
    const int k, const double alpha, const double* const a, const int lda,
    const double* const b, const int ldb, const double beta, double* const c,
    const int ldc);

#ifdef HAVE_MAGMA
template void LAU_D::MPgemm<float, float, float>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const float* const a, const int lda,
    const float* const b, const int ldb, const double beta, float* const c,
    const int ldc);
template void LAU_D::MPgemm<double, float, double>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const double* const a, const int lda,
    const float* const b, const int ldb, const double beta, double* const c,
    const int ldc);
template void LAU_D::MPgemm<float, double, float>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const float* const a, const int lda,
    const double* const b, const int ldb, const double beta, float* const c,
    const int ldc);
template void LAU_D::MPgemm<double, double, float>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const double* const a, const int lda,
    const double* const b, const int ldb, const double beta, float* const c,
    const int ldc);
template void LAU_D::MPgemm<float, float, double>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const float* const a, const int lda,
    const float* const b, const int ldb, const double beta, double* const c,
    const int ldc);
template void LAU_D::MPgemmNN<float, double, float>(const int m, const int n,
    const int k, const double alpha, const float* const a, const int lda,
    const double* const b, const int ldb, const double beta, float* const c,
    const int ldc);
template void LAU_D::MPgemmNN<double, double, double>(const int m, const int n,
    const int k, const double alpha, const double* const a, const int lda,
    const double* const b, const int ldb, const double beta, double* const c,
    const int ldc);
template void LAU_D::MPscal(const int len, const double scal, double* ptr);
template void LAU_D::MPscal(const int len, const double scal, float* ptr);
template double LAU_D::MPdot<double, double>(
    const int len, const double* const xptr, const double* const yptr);
template double LAU_D::MPdot<float, float>(const int len,
    const float* __restrict__ xptr, const float* __restrict__ yptr);
template double LAU_D::MPdot<double, float>(
    const int len, const double* const xptr, const float* const yptr);
template double LAU_D::MPdot<float, double>(
    const int len, const float* const xptr, const double* const yptr);
template void LAU_D::MPaxpy<float, double>(const int len, const double scal,
    const float* __restrict__ xptr, double* __restrict__ yptr);
template void LAU_D::MPaxpy<float, float>(const int len, const double scal,
    const float* __restrict__ xptr, float* __restrict__ yptr);
#endif
