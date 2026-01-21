// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_MPUTILS_H
#define MGMOL_MPUTILS_H

#include "MGmol_blas1.h"
#include "memory_space.h"

/* scal */
/* scalar times vector for use in templates.
 * Calls blas dscal or sscal depending on arguments
 */
void Tscal(const int len, const double scal, double* dptr);
void Tscal(const int len, const float scal, float* dptr);

/* nrm2 */
/* Vector 2-norm for use in templates.
 * Calls blas dnrm2 or snrm2 depending on arguments
 */
double Tnrm2(const int len, const double* const dptr);
float Tnrm2(const int len, const float* const dptr);
/* mixed-precision vector 2-norm. Accumulates results
 * in double precision and stores as double precision.
 */
double MPnrm2(const int len, const float* const dptr);

/* dot */
/* Vector dot-product for use in templates.
 * Calls blas ddot or sdot depending on arguments
 */
double Tdot(const int len, const double* const xptr, const double* const yptr);
float Tdot(const int len, const float* const xptr, const float* const yptr);
/* axpy */
/* Scalar times vector plus vector (AXPY) for use in templates.
 * Calls blas daxpy or saxpy depending on arguments
 */
void Taxpy(const int len, double scal, const double* const xptr, double* yptr);
void Taxpy(const int len, float scal, const float* const xptr, float* yptr);

/* copy */
/* Vector copy for use in templates.
 * Calls blas scopy or dcopy depending on arguments.
 */
inline void Tcopy(const int* const len, const double* const x,
    const int* const incx, double* y, const int* const incy)
{
    MemorySpace::assert_is_host_ptr(x);
    MemorySpace::assert_is_host_ptr(x);
    DCOPY(len, x, incx, y, incy);
}
inline void Tcopy(const int* const len, const float* const x,
    const int* const incx, float* y, const int* const incy)
{
    MemorySpace::assert_is_host_ptr(x);
    MemorySpace::assert_is_host_ptr(x);
    SCOPY(len, x, incx, y, incy);
}

/* syrk */
/* Symmetric rank-k update: C = alpha*A**T*A + beta*C  or C = alpha*A*A**T +
 * beta*C Calls blas ssyrk or dsyrk depending on arguments
 */
void Tsyrk(const char uplo, const char trans, const int n, const int k,
    const double alpha, const double* const a, const int lda, const double beta,
    const double* c, const int ldc);
void Tsyrk(const char uplo, const char trans, const int n, const int k,
    const float alpha, const float* const a, const int lda, const float beta,
    float* c, const int ldc);

void Tgemv(const char trans, const int m, const int n, const double alpha,
    const double* const a, const int lda, const double* const x, const int incx,
    const double beta, double* const y, const int incy);
void Tgemv(const char trans, const int m, const int n, const float alpha,
    const float* const a, const int lda, const float* const x, const int incx,
    const float beta, float* const y, const int incy);

int Tiamax(const int*, const double* const a, const int*);
int Tiamax(const int*, const float* const a, const int*);

template <typename MemorySpaceType>
struct LinearAlgebraUtils
{
    template <typename T1, typename T2, typename T3>
    static void MPgemm(const char transa, const char transb, const int m,
        const int n, const int k, const double alpha, const T1* const a,
        const int lda, const T2* const b, const int ldb, const double beta,
        T3* const c, const int ldc);

    template <typename T1, typename T2, typename T3>
    static void MPgemmNN(const int m, const int n, const int k,
        const double alpha, const T1* const a, const int lda, const T2* const b,
        const int ldb, const double beta, T3* const c, const int ldc);

    template <typename T1, typename T2, typename T3>
    static void MPgemmTN(const int m, const int n, const int k,
        const double alpha, const T1* const a, const int lda, const T2* const b,
        const int ldb, const double beta, T3* const c, const int ldc);

    /* mixed-precision scalar times vector. Accumulates results
     * in double precision and stores as single precision.
     */
    static void MPscal(const int len, const double scal, double* dptr);
    static void MPscal(const int len, const double scal, float* dptr);

    /* mixed-precision vector dot-product. Accumulates results
     * in double precision and stores as double precision.
     */
    template <typename T1, typename T2>
    static double MPdot(
        const int len, const T1* const xptr, const T2* const yptr);

    /* mixed-precision vector times scalar plus vector. Accumulates results
     * in double precision and stores in single precision.
     */
    template <typename T0, typename T1, typename T2>
    static void MPaxpy(const int len, T0 scal, const T1* __restrict__ xptr,
        T2* __restrict__ yptr);

    static void MPsyrk(const char uplo, const char trans, const int n,
        const int k, const double alpha, const double* const a, const int lda,
        const double beta, double* c, const int ldc);
    template <typename T1, typename T2>
    static void MPsyrk(const char uplo, const char trans, const int n,
        const int k, const double alpha, const T1* const a, const int lda,
        const double beta, T2* c, const int ldc);
};

/* trsm */
void Ttrsm(const char, const char, const char, const char, const int, const int,
    const double, const double* const, const int, double* const, const int);
void Ttrsm(const char, const char, const char, const char, const int, const int,
    const float, const float* const, const int, float* const, const int);

void MPcpy(double* const dest, const double* const src, const int n);
void MPcpy(float* const dest, const float* const src, const int n);
void MPcpy(
    double* __restrict__ dest, const float* __restrict__ src, const int n);
void MPcpy(
    float* __restrict__ dest, const double* __restrict__ src, const int n);

#endif
