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

#include "MGmol_blas1.h"
#include "gemm_impl.h"
#include "syrk_impl.h"

#ifdef MGMOL_USE_BLIS
#include <blis.h>
#else
#include "blas2_c.h"
#include "blas3_c.h"
#endif

#include <cassert>
#include <iostream>
#include <vector>

Timer dgemm_tm("dgemm");
Timer sgemm_tm("sgemm");
Timer bligemm_tm("bligemm");

Timer dsyrk_tm("dsyrk");
Timer ssyrk_tm("ssyrk");

Timer ttdot_tm("ttdot");

// Timers for hand written loops
Timer loopdot_tm("loopdot");
Timer loopaxpy_tm("loopaxpy");
Timer loopscal_tm("loopscal");
Timer loopcp_tm("loopcp");

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
    MemorySpace::assert_is_host_ptr(dptr);

    const int one = 1;
    DSCAL(&len, &scal, dptr, &one);
}

template <>
void LAU_H::MPscal(const int len, const double scal, float* dptr)
{
    loopscal_tm.start();

    MemorySpace::assert_is_host_ptr(dptr);

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

    loopscal_tm.stop();
}

// MemorySpace::Device
#ifdef HAVE_MAGMA
template <>
void LAU_D::MPscal(const int len, const double scal, double* dptr)
{
    MemorySpace::assert_is_dev_ptr(dptr);

    int const increment   = 1;
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();
    magma_dscal(len, scal, dptr, increment, magma_singleton.queue_);
}

template <>
void LAU_D::MPscal(const int len, const double scal, float* dptr)
{
    MemorySpace::assert_is_dev_ptr(dptr);

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
        MemorySpace::copy_to_host(dptr, len, dptr_alias.get());
#endif

        MGMOL_PARALLEL_FOR(dptr_alias)
        for (int k = 0; k < len; k++)
        {
            double val    = static_cast<double>(dptr_alias[k]);
            dptr_alias[k] = static_cast<float>(scal * val);
        }

#ifndef HAVE_OPENMP_OFFLOAD
        MemorySpace::copy_to_dev(dptr_alias.get(), len, dptr);
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
    MemorySpace::assert_is_host_ptr(xptr);
    MemorySpace::assert_is_host_ptr(yptr);

    const int one = 1;
    return DDOT(&len, xptr, &one, yptr, &one);
}

template <>
template <typename T1, typename T2>
double LAU_H::MPdot(
    const int len, const T1* __restrict__ xptr, const T2* __restrict__ yptr)
{
    MemorySpace::assert_is_host_ptr(xptr);
    MemorySpace::assert_is_host_ptr(yptr);

    loopdot_tm.start();

    double dot = 0.;
    for (int k = 0; k < len; k++)
    {
        double val1 = static_cast<double>(xptr[k]);
        double val2 = static_cast<double>(yptr[k]);
        dot += val1 * val2;
    }

    loopdot_tm.stop();

    return dot;
}

// MemorySpace::Device
#ifdef HAVE_MAGMA
template <>
template <>
double LAU_D::MPdot<double, double>(
    const int len, const double* const xptr, const double* const yptr)
{
    MemorySpace::assert_is_dev_ptr(xptr);
    MemorySpace::assert_is_dev_ptr(yptr);

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
    MemorySpace::assert_is_dev_ptr(xptr);
    MemorySpace::assert_is_dev_ptr(yptr);

#ifndef HAVE_OPENMP_OFFLOAD
    std::unique_ptr<T1[], void (*)(T1*)> xptr_host(
        MemoryH<T1>::allocate(len), MemoryH<T1>::free);
    std::unique_ptr<T2[], void (*)(T2*)> yptr_host(
        MemoryH<T2>::allocate(len), MemoryH<T2>::free);
    MemorySpace::copy_to_host(xptr, len, xptr_host.get());
    MemorySpace::copy_to_host(yptr, len, yptr_host.get());
    return LAU_H::MPdot(len, xptr_host.get(), yptr_host.get());
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
template <>
void LAU_H::MPaxpy(const int len, double scal, const double* __restrict__ xptr,
    double* __restrict__ yptr)
{
    MemorySpace::assert_is_host_ptr(xptr);
    MemorySpace::assert_is_host_ptr(yptr);

    const int one = 1;
    DAXPY(&len, &scal, xptr, &one, yptr, &one);
}

template <>
template <>
void LAU_H::MPaxpy(const int len, float scal, const float* __restrict__ xptr,
    float* __restrict__ yptr)
{
    MemorySpace::assert_is_host_ptr(xptr);
    MemorySpace::assert_is_host_ptr(yptr);

    const int one = 1;
    SAXPY(&len, &scal, xptr, &one, yptr, &one);
}

template <>
template <typename T0, typename T1, typename T2>
void LAU_H::MPaxpy(
    const int len, T0 scal, const T1* __restrict__ xptr, T2* __restrict__ yptr)
{
    loopaxpy_tm.start();

    MemorySpace::assert_is_host_ptr(xptr);
    MemorySpace::assert_is_host_ptr(yptr);
#pragma omp parallel for simd
    for (int k = 0; k < len; k++)
    {
        yptr[k] += static_cast<T2>(scal * static_cast<double>(xptr[k]));
    }

    loopaxpy_tm.stop();
}

// MemorySpace::Device
#ifdef HAVE_MAGMA
template <>
void LAU_D::MPaxpy(const int len, double scal, const double* __restrict__ xptr,
    double* __restrict__ yptr)
{
    MemorySpace::assert_is_dev_ptr(xptr);
    MemorySpace::assert_is_dev_ptr(yptr);

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
    MemorySpace::assert_is_dev_ptr(xptr);
    MemorySpace::assert_is_dev_ptr(yptr);

#ifndef HAVE_OPENMP_OFFLOAD
    std::unique_ptr<T1[], void (*)(T1*)> xptr_host(
        MemoryH<T1>::allocate(len), MemoryH<T1>::free);
    std::unique_ptr<T2[], void (*)(T2*)> yptr_host(
        MemoryH<T2>::allocate(len), MemoryH<T2>::free);
    MemorySpace::copy_to_host(xptr, len, xptr_host.get());
    MemorySpace::copy_to_host(yptr, len, yptr_host.get());
    LAU_H::MPaxpy(len, scal, xptr_host.get(), yptr_host.get());
    MemorySpace::copy_to_dev(yptr_host.get(), len, yptr);
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

///////////////////////////////
////          MPsyrk         //
///////////////////////////////
// MemorySpace::Host
template <>
void LAU_H::MPsyrk(const char uplo, const char trans, const int n, const int k,
    const double alpha, const double* const a, const int lda, const double beta,
    double* c, const int ldc)
{
    MemorySpace::assert_is_host_ptr(a);
    MemorySpace::assert_is_host_ptr(c);

    dsyrk_tm.start();
    DSYRK(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
    dsyrk_tm.stop();
}

template <>
template <typename T1, typename T2>
void LAU_H::MPsyrk(const char uplo, const char trans, const int n, const int k,
    const double alpha, const T1* const a, const int lda, const double beta,
    T2* c, const int ldc)
{
    MemorySpace::assert_is_host_ptr(a);
    MemorySpace::assert_is_host_ptr(c);

    syrk_impl(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}

// MemorySpace::Device
#ifdef HAVE_MAGMA
template <>
void LAU_D::MPsyrk(const char uplo, const char trans, const int n, const int k,
    const double alpha, const double* const a, const int lda, const double beta,
    double* c, const int ldc)
{
    MemorySpace::assert_is_dev_ptr(a);
    MemorySpace::assert_is_dev_ptr(c);

    dsyrk_tm.start();
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();
    magma_dsyrk(magma_uplo_const(uplo), magma_trans_const(trans), n, k, alpha,
        a, lda, beta, c, ldc, magma_singleton.queue_);
    dsyrk_tm.stop();
}

template <>
void LAU_D::MPsyrk(const char uplo, const char trans, const int n, const int k,
    const double alpha, const float* const a, const int lda, const double beta,
    float* c, const int ldc)
{
    MemorySpace::assert_is_dev_ptr(a);
    MemorySpace::assert_is_dev_ptr(c);

    // Move the data on the host
    const int size_a = lda * k;
    std::vector<float> a_host(size_a);
    MemorySpace::copy_to_host(a, a_host);

    const int size_c = ldc * n;
    std::vector<float> c_host(size_c);
    MemorySpace::copy_to_host(c, c_host);

    // Call the host function
    LAU_H::MPsyrk(
        uplo, trans, n, k, alpha, a_host.data(), lda, beta, c_host.data(), ldc);

    // Move the result back to the
    MemorySpace::copy_to_dev(c_host, c);
}

template <>
template <typename T1, typename T2>
void LAU_D::MPsyrk(const char uplo, const char trans, const int n, const int k,
    const double alpha, const T1* const a, const int lda, const double beta,
    T2* c, const int ldc)
{
    MemorySpace::assert_is_dev_ptr(a);
    MemorySpace::assert_is_dev_ptr(c);

    // Move the data on the host
    const int size_a = lda * k;
    std::vector<T1> a_host(size_a);
    MemorySpace::copy_to_host(a, a_host);

    const int size_c = ldc * n;
    std::vector<T2> c_host(size_c);
    MemorySpace::copy_to_host(c, c_host);

    // Call the host function
    LAU_H::MPsyrk(
        uplo, trans, n, k, alpha, a_host.data(), lda, beta, c_host.data(), ldc);

    // Move the result back to the
    MemorySpace::copy_to_dev(c_host, c);
}
#endif

void Tscal(const int len, const double scal, double* dptr)
{
    MemorySpace::assert_is_host_ptr(dptr);
    const int one = 1;
    DSCAL(&len, &scal, dptr, &one);
}
void Tscal(const int len, const float scal, float* dptr)
{
    MemorySpace::assert_is_host_ptr(dptr);
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

int Tiamax(const int* n, const float* const a, const int* incx)
{
    return ISAMAX(n, a, incx);
}

int Tiamax(const int* n, const double* const a, const int* incx)
{
    return IDAMAX(n, a, incx);
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
    MemorySpace::assert_is_host_ptr(a);
    MemorySpace::assert_is_host_ptr(b);
    MemorySpace::assert_is_host_ptr(c);

    gemm_impl(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

// input/output in double, computation in double
template <>
template <>
void LAU_H::MPgemm<double, double, double>(const char transa, const char transb,
    const int m, const int n, const int k, const double alpha,
    const double* const a, const int lda, const double* const b, const int ldb,
    const double beta, double* const c, const int ldc)
{
    MemorySpace::assert_is_host_ptr(a);
    MemorySpace::assert_is_host_ptr(b);
    MemorySpace::assert_is_host_ptr(c);

    dgemm_tm.start();
    DGEMM(
        &transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    dgemm_tm.stop();
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
    MemorySpace::assert_is_dev_ptr(a);
    MemorySpace::assert_is_dev_ptr(b);
    MemorySpace::assert_is_dev_ptr(c);

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
    MemorySpace::assert_is_dev_ptr(a);
    MemorySpace::assert_is_dev_ptr(b);
    MemorySpace::assert_is_dev_ptr(c);

    dgemm_tm.start();
    // Transform char to magma_trans_t
    auto convert_to_magma_trans = [](const char trans)
    {
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

template <>
template <>
void LAU_H::MPgemmNN(const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const double* const b, const int ldb,
    const double beta, float* const c, const int ldc)
{
    MemorySpace::assert_is_host_ptr(a);
    MemorySpace::assert_is_host_ptr(b);
    MemorySpace::assert_is_host_ptr(c);

#ifdef MGMOL_USE_BLIS
    bligemm_tm.start();

    // Create matrix objects
    // When storing by columns, the row stride is 1
    // When storing by columns, the column stride is also sometimes called the
    // leading dimension
    obj_t A;
    bli_obj_create_with_attached_buffer(
        BLIS_FLOAT, m, k, const_cast<float*>(a), 1, lda, &A);

    obj_t B;
    bli_obj_create_with_attached_buffer(
        BLIS_DOUBLE, k, n, const_cast<double*>(b), 1, ldb, &B);

    obj_t C;
    bli_obj_create_with_attached_buffer(
        BLIS_FLOAT, m, n, const_cast<float*>(c), 1, ldc, &C);

    obj_t bli_alpha;
    bli_obj_create_1x1(BLIS_DOUBLE, &bli_alpha);
    bli_setsc(alpha, 0., &bli_alpha);

    obj_t bli_beta;
    bli_obj_create_1x1(BLIS_DOUBLE, &bli_beta);
    bli_setsc(beta, 0., &bli_beta);

    // accumulate results in double precision
    bli_obj_set_comp_prec(BLIS_DOUBLE_PREC, &C);

    bli_gemm(&bli_alpha, &A, &B, &bli_beta, &C);

    // Clean up BLIS objects
    bli_obj_free(&bli_alpha);
    bli_obj_free(&bli_beta);

    bligemm_tm.stop();
#else
    char transa = 'n';
    char transb = 'n';

    LAU_H::MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#endif
}

// input in float, computation in double
template <>
template <>
void LAU_H::MPgemmTN(const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const float* const b, const int ldb,
    const double beta, double* const c, const int ldc)
{
    // std::cout << "LAU_H::MPgemmTN" << std::endl;
    MemorySpace::assert_is_host_ptr(a);
    MemorySpace::assert_is_host_ptr(b);
    MemorySpace::assert_is_host_ptr(c);

#ifdef MGMOL_USE_BLIS
    bligemm_tm.start();

    // Create matrix objects
    // When storing by columns, the row stride is 1
    // When storing by columns, the column stride is also sometimes called the
    // leading dimension
    obj_t A;
    bli_obj_create_with_attached_buffer(
        BLIS_FLOAT, k, m, const_cast<float*>(a), 1, lda, &A);
    bli_obj_toggle_trans(&A);

    obj_t B;
    bli_obj_create_with_attached_buffer(
        BLIS_FLOAT, k, n, const_cast<float*>(b), 1, ldb, &B);
    obj_t C;
    bli_obj_create_with_attached_buffer(
        BLIS_DOUBLE, m, n, const_cast<double*>(c), 1, ldc, &C);

    obj_t bli_alpha;
    bli_obj_create_1x1(BLIS_DOUBLE, &bli_alpha);
    bli_setsc(alpha, 0., &bli_alpha);

    obj_t bli_beta;
    bli_obj_create_1x1(BLIS_DOUBLE, &bli_beta);
    bli_setsc(beta, 0., &bli_beta);

    // accumulate results in double precision
    // dafault: precision of C
    bli_obj_set_comp_prec(BLIS_DOUBLE_PREC, &C);
    bli_gemm(&bli_alpha, &A, &B, &bli_beta, &C);

    // Clean up BLIS objects
    bli_obj_free(&bli_alpha);
    bli_obj_free(&bli_beta);

    bligemm_tm.stop();
#else
    char transa = 't';
    char transb = 'n';

    LAU_H::MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#endif
}

// input in float, computation in double
template <>
template <>
void LAU_H::MPgemmTN(const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const float* const b, const int ldb,
    const double beta, float* const c, const int ldc)
{
    char transa = 't';
    char transb = 'n';

    LAU_H::MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

template <typename MemorySpaceType>
template <typename T1, typename T2, typename T3>
void LinearAlgebraUtils<MemorySpaceType>::MPgemmTN(const int m, const int n,
    const int k, const double alpha, const T1* const a, const int lda,
    const T2* const b, const int ldb, const double beta, T3* const c,
    const int ldc)
{
    // if(onpe0)cout<<"template MPgemmNN..."<<endl;
    char transa = 't';
    char transb = 'n';
    LAU_H::MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
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
    loopcp_tm.start();

    for (int i = 0; i < n; i++)
        dest[i] = src[i];

    loopcp_tm.stop();
}
void MPcpy(
    float* __restrict__ dest, const double* __restrict__ src, const int n)
{
    loopcp_tm.start();

    for (int i = 0; i < n; i++)
        dest[i] = src[i];

    loopcp_tm.stop();
}

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
template void LAU_H::MPgemm<float, float, float>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const float* const a, const int lda,
    const float* const b, const int ldb, const double beta, float* const c,
    const int ldc);

// template void LAU_H::MPgemmNN<float, double, float>(const int m, const int n,
//    const int k, const double alpha, const float* const a, const int lda,
//    const double* const b, const int ldb, const double beta, float* const c,
//    const int ldc);
template void LAU_H::MPgemmNN<double, double, double>(const int m, const int n,
    const int k, const double alpha, const double* const a, const int lda,
    const double* const b, const int ldb, const double beta, double* const c,
    const int ldc);
template double LAU_H::MPdot<float, float>(const int len,
    const float* __restrict__ xptr, const float* __restrict__ yptr);
template double LAU_H::MPdot<double, float>(
    const int len, const double* const xptr, const float* const yptr);
template double LAU_H::MPdot<float, double>(
    const int len, const float* const xptr, const double* const yptr);
template void LAU_H::MPaxpy<double, float, double>(const int len,
    const double scal, const float* __restrict__ xptr,
    double* __restrict__ yptr);
template void LAU_H::MPaxpy<float, float, double>(const int len,
    const float scal, const float* __restrict__ xptr,
    double* __restrict__ yptr);
template void LAU_H::MPaxpy<double, float, float>(const int len,
    const double scal, const float* __restrict__ xptr,
    float* __restrict__ yptr);

template void LAU_H::MPsyrk<double, float>(const char uplo, const char trans,
    const int n, const int k, const double alpha, const double* const a,
    const int lda, const double beta, float* c, const int ldc);
template void LAU_H::MPsyrk<float, double>(const char uplo, const char trans,
    const int n, const int k, const double alpha, const float* const a,
    const int lda, const double beta, double* c, const int ldc);
template void LAU_H::MPsyrk<float, float>(const char uplo, const char trans,
    const int n, const int k, const double alpha, const float* const a,
    const int lda, const double beta, float* c, const int ldc);

template void LAU_H::MPgemmTN<double, double, double>(const int m, const int n,
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
template void LAU_D::MPsyrk<double, float>(const char uplo, const char trans,
    const int n, const int k, const double alpha, const double* const a,
    const int lda, const double beta, float* c, const int ldc);
template void LAU_D::MPsyrk<float, double>(const char uplo, const char trans,
    const int n, const int k, const double alpha, const float* const a,
    const int lda, const double beta, double* c, const int ldc);
#endif
