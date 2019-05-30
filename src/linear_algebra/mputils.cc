// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "mputils.h"
#include "MPIdata.h"
#include "Timer.h"

#include <iostream>
using namespace std;


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
void MPscal(const int len, const double scal, double* dptr)
{
    const int one = 1;
    DSCAL(&len, &scal, dptr, &one);
}

void MPscal(const int len, const double scal, float* dptr)
{
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
            double val = (double)dptr[k];
            dptr[k]    = (float)(scal * val);
        }
    }
}

double Tnrm2(const int len, const double* const dptr)
{
    const int one = 1;
    double nrm;
    nrm = dnrm2(&len, dptr, &one);

    return nrm;
}
float Tnrm2(const int len, const float* const dptr)
{
    const int one = 1;
    float nrm;
    nrm = snrm2(&len, dptr, &one);

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
    dot = ddot(&len, xptr, &one, yptr, &one);

    return dot;
}
float Tdot(const int len, const float* const xptr, const float* const yptr)
{
    const int one = 1;
    float dot;
    dot = sdot(&len, xptr, &one, yptr, &one);

    return dot;
}
double MPdot(const int len, const float* __restrict__ xptr,
    const float* __restrict__ yptr)
{
    mpdot_tm.start();

#ifdef BGQ
    const int one = 1;
    double dot    = (double)sdot(&len, xptr, &one, yptr, &one);
#else
    double dot = 0.;
    for (int k = 0; k < len; k++)
    {
        double val1 = (double)xptr[k];
        double val2 = (double)yptr[k];
        dot += val1 * val2;
    }
#endif

    mpdot_tm.stop();

    return dot;
}
double MPdot(const int len, const double* const xptr, const double* const yptr)
{
    const int one = 1;
    return ddot(&len, xptr, &one, yptr, &one);
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
void MPaxpy(const int len, double scal, const double* __restrict__ xptr,
    double* __restrict__ yptr)
{
    const int one = 1;
    DAXPY(&len, &scal, xptr, &one, yptr, &one);
}

void Ttrsm(const char side, const char uplo, const char transa, const char diag,
    const int m, const int n, const double alpha, const double* const a,
    const int lda, double* const b, const int ldb)
{
    dtrsm(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

void Ttrsm(const char side, const char uplo, const char transa, const char diag,
    const int m, const int n, const float alpha, const float* const a,
    const int lda, float* const b, const int ldb)
{
    strsm(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

void Tsyrk(const char uplo, const char trans, const int n, const int k,
    const double alpha, const double* const a, const int lda, const double beta,
    double* c, const int ldc)
{
    dsyrk_tm.start();
    dsyrk(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
    dsyrk_tm.stop();
}
void Tsyrk(const char uplo, const char trans, const int n, const int k,
    const float alpha, const float* const a, const int lda, const float beta,
    float* c, const int ldc)
{
    ssyrk_tm.start();
    ssyrk(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
    ssyrk_tm.stop();
}

void MPsyrk(const char uplo, const char trans, const int n, const int k,
    const double alpha, const double* const a, const int lda, const double beta,
    double* c, const int ldc)
{
    dsyrk_tm.start();
    dsyrk(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
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
        double* buff = new double[n];
        if (uplo == 'U' || uplo == 'u')
        {
            for (int j = 0; j < n; j++)
            {
                const int len = j + 1;
                memset(buff, 0, len * sizeof(double));
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const float* colL = a + lda * l;
                    /* get multiplier */
                    double mult
                        = (double)(alpha
                                   * colL[j]); // same as alpha * a[lda*l + j];
                    MPaxpy(len, mult, colL, buff);
                }
                /* Update col j of upper part of matrix C. */
                /* Get pointer to beginning of column j in C. */
                float* cj = c + ldc * j;
                MPscal(len, beta, cj);
                for (int i = 0; i < len; i++)
                    cj[i] += (float)buff[i];
            }
        }
        else /* uplo = 'L' or 'l' */
        {
            for (int j = 0; j < n; j++)
            {
                const int len = n - (j + 1);
                memset(buff, 0, len * sizeof(double));
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const float* colL = a + lda * l + j;
                    /* get multiplier */
                    double mult
                        = (double)(alpha
                                   * colL[0]); // same as alpha * a[lda*l + j];
                    MPaxpy(len, mult, colL, buff);
                }
                /* Update col j of upper part of matrix C. */
                /* Get pointer to beginning of column j in C. */
                float* cj = c + ldc * j + j;
                MPscal(len, beta, cj);
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
                    double bc                    = (double)c[pos] * beta;
                    c[pos] = (float)(alpha * MPdot(k, ai, aj) + bc);
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
                    double bc                    = (double)c[pos] * beta;
                    c[pos] = (float)(alpha * MPdot(k, ai, aj) + bc);
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

void MPgemm(const char transa, const char transb, const int m, const int n,
    const int k, const double alpha, const double* const a, const int lda,
    const double* const b, const int ldb, const double beta, double* const c,
    const int ldc)
{
    dgemm_tm.start();
    DGEMM(
        &transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    dgemm_tm.stop();
}

// input/output in float, computation in double
void MPgemmNN(const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const double* const b, const int ldb,
    const double beta, float* const c, const int ldc)
{
    // if(onpe0)cout<<"MPgemmNN..."<<endl;

    char transa = 'n';
    char transb = 'n';

    MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

// input in float, computation in double
void MPgemmTN(const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const float* const b, const int ldb,
    const double beta, double* const c, const int ldc)
{
    char transa = 't';
    char transb = 'n';

    MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

// input in float, computation in double
void MPgemmTN(const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const float* const b, const int ldb,
    const double beta, float* const c, const int ldc)
{
    char transa = 't';
    char transb = 'n';

    MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

// input/output in float, computation in double
void MPgemm(const char transa, const char transb, const int m, const int n,
    const int k, const double alpha, const float* const a, const int lda,
    const float* const b, const int ldb, const double beta, float* const c,
    const int ldc)
{
    mpgemm_tm.start();

    if (beta == 1. && (alpha == 0. || m == 0 || n == 0 || k == 0)) return;

    /* case transb == 'N' and transa == 'N' */
    if (transb == 'N' || transb == 'n')
    {
        if (transa == 'N' || transa == 'n')
        {
            /* buffer to hold accumulation in double */
            double* buff = new double[m];
            for (int j = 0; j < n; j++)
            {
                memset(buff, 0, m * sizeof(double));
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const float* colL = a + lda * l;
                    /* get multiplier */
                    double mult = (double)(alpha * b[ldb * j + l]);
                    MPaxpy(m, mult, colL, buff);
                }
                /* Update col j of of result matrix C. */
                /* Get pointer to beginning of column j in C. */
                float* cj = c + ldc * j;
                MPscal(m, beta, cj);
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
            double* buff = new double[m];
            for (int j = 0; j < n; j++)
            {
                memset(buff, 0, m * sizeof(double));
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const float* colL = a + lda * l;
                    /* get multiplier */
                    double mult = (double)(alpha * b[ldb * l + j]);
                    MPaxpy(m, mult, colL, buff);
                }
                /* Update col j of of result matrix C. */
                /* Get pointer to beginning of column j in C. */
                float* cj = c + ldc * j;
                MPscal(m, beta, cj);
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

/////// additional calls ... may be removed later if unused
template <typename T1, typename T2>
double MPdot(
    const int len, const T1* __restrict__ xptr, const T2* __restrict__ yptr)
{
    ttdot_tm.start();

    double dot = 0.;
    for (int k = 0; k < len; k++)
    {
        double val1 = (double)xptr[k];
        double val2 = (double)yptr[k];
        dot += val1 * val2;
    }

    ttdot_tm.stop();

    return dot;
}
template <typename T1, typename T2>
void MPaxpy(const int len, double scal, const T1* __restrict__ xptr,
    T2* __restrict__ yptr)
{
    for (int k = 0; k < len; k++)
    {
        yptr[k] += (T2)(scal * (double)xptr[k]);
    }
}

template <typename T1, typename T2, typename T3>
void MPgemmNN(const int m, const int n, const int k, const double alpha,
    const T1* const a, const int lda, const T2* const b, const int ldb,
    const double beta, T3* const c, const int ldc)
{
    // if(onpe0)cout<<"template MPgemmNN..."<<endl;
    char transa = 'n';
    char transb = 'n';
    MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

template <typename T1, typename T2, typename T3>
void MPgemmTN(const int m, const int n, const int k, const double alpha,
    const T1* const a, const int lda, const T2* const b, const int ldb,
    const double beta, T3* const c, const int ldc)
{
    // if(onpe0)cout<<"template MPgemmNN..."<<endl;
    char transa = 't';
    char transb = 'n';
    MPgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

template <typename T1, typename T2, typename T3>
void MPgemm(const char transa, const char transb, const int m, const int n,
    const int k, const double alpha, const T1* const a, const int lda,
    const T2* const b, const int ldb, const double beta, T3* const c,
    const int ldc)
{

    tttgemm_tm.start();
    // if(onpe0)cout<<"template MPgemm..."<<endl;

    if (beta == 1. && (alpha == 0. || m == 0 || n == 0 || k == 0)) return;

    /* case transb == 'N' and transa == 'N' */
    if (transb == 'N' || transb == 'n')
    {
        if (transa == 'N' || transa == 'n')
        {
            /* buffer to hold accumulation in double */
            double* buff = new double[m];
            for (int j = 0; j < n; j++)
            {
                memset(buff, 0, m * sizeof(double));
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const T1* colL = a + lda * l;
                    /* get multiplier */
                    double mult = (double)(alpha * b[ldb * j + l]);
                    MPaxpy(m, mult, colL, buff);
                }
                /* Update col j of of result matrix C. */
                /* Get pointer to beginning of column j in C. */
                T3* cj = c + ldc * j;
                MPscal(m, beta, cj);
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
                    const int pos             = ldc * j + i;
                    double bc                 = (double)c[pos] * beta;
                    const T1* __restrict__ ai = a + lda * i;
                    c[pos] = (T3)(alpha * MPdot(k, ai, bj) + bc);
                }
            }
        }
    }
    else /* transb == 'T'/'C' */
    {
        if (transa == 'N' || transa == 'n')
        {
            /* buffer to hold accumulation in double */
            double* buff = new double[m];
            for (int j = 0; j < n; j++)
            {
                memset(buff, 0, m * sizeof(double));
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const T1* colL = a + lda * l;
                    /* get multiplier */
                    double mult = (double)(alpha * b[ldb * l + j]);
                    MPaxpy(m, mult, colL, buff);
                }
                /* Update col j of of result matrix C. */
                /* Get pointer to beginning of column j in C. */
                T3* cj = c + ldc * j;
                MPscal(m, beta, cj);
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
        double* buff = new double[n];
        if (uplo == 'U' || uplo == 'u')
        {
            for (int j = 0; j < n; j++)
            {
                const int len = j + 1;
                memset(buff, 0, len * sizeof(double));
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const T1* colL = a + lda * l;
                    /* get multiplier */
                    double mult
                        = (double)(alpha
                                   * colL[j]); // same as alpha * a[lda*l + j];
                    MPaxpy(len, mult, colL, buff);
                }
                /* Update col j of upper part of matrix C. */
                /* Get pointer to beginning of column j in C. */
                T2* cj = c + ldc * j;
                MPscal(len, beta, cj);
                for (int i = 0; i < len; i++)
                    cj[i] += (T2)buff[i];
            }
        }
        else /* uplo = 'L' or 'l' */
        {
            for (int j = 0; j < n; j++)
            {
                const int len = n - (j + 1);
                memset(buff, 0, len * sizeof(double));
                for (int l = 0; l < k; l++)
                {
                    /* pointer to beginning of column l in matrix a */
                    const T1* colL = a + lda * l + j;
                    /* get multiplier */
                    double mult
                        = (double)(alpha
                                   * colL[0]); // same as alpha * a[lda*l + j];
                    MPaxpy(len, mult, colL, buff);
                }
                /* Update col j of upper part of matrix C. */
                /* Get pointer to beginning of column j in C. */
                T2* cj = c + ldc * j + j;
                MPscal(len, beta, cj);
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
                    double bc                 = (double)c[pos] * beta;
                    c[pos] = (T2)(alpha * MPdot(k, ai, aj) + bc);
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
                    double bc                 = (double)c[pos] * beta;
                    c[pos] = (T2)(alpha * MPdot(k, ai, aj) + bc);
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

template double MPdot<double, float>(const int len,
    const double* __restrict__ xptr, const float* __restrict__ yptr);
template double MPdot<float, double>(const int len,
    const float* __restrict__ xptr, const double* __restrict__ yptr);

template void MPaxpy<float, double>(const int len, const double scal,
    const float* __restrict__ xptr, double* __restrict__ yptr);
template void MPaxpy<float, float>(const int len, const double scal,
    const float* __restrict__ xptr, float* __restrict__ yptr);

template void MPsyrk<double, float>(const char uplo, const char trans,
    const int n, const int k, const double alpha, const double* const a,
    const int lda, const double beta, float* c, const int ldc);
template void MPsyrk<float, double>(const char uplo, const char trans,
    const int n, const int k, const double alpha, const float* const a,
    const int lda, const double beta, double* c, const int ldc);

template void MPgemm<double, float, double>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const double* const a, const int lda,
    const float* const b, const int ldb, const double beta, double* const c,
    const int ldc);
template void MPgemm<float, double, float>(const char transa, const char transb,
    const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const double* const b, const int ldb,
    const double beta, float* const c, const int ldc);
template void MPgemm<double, double, float>(const char transa,
    const char transb, const int m, const int n, const int k,
    const double alpha, const double* const a, const int lda,
    const double* const b, const int ldb, const double beta, float* const c,
    const int ldc);
template void MPgemm<float, float, double>(const char transa, const char transb,
    const int m, const int n, const int k, const double alpha,
    const float* const a, const int lda, const float* const b, const int ldb,
    const double beta, double* const c, const int ldc);
template void MPgemmNN<double, double, double>(const int m, const int n,
    const int k, const double alpha, const double* const a, const int lda,
    const double* const b, const int ldb, const double beta, double* const c,
    const int ldc);
template void MPgemmTN<double, double, double>(const int m, const int n,
    const int k, const double alpha, const double* const a, const int lda,
    const double* const b, const int ldb, const double beta, double* const c,
    const int ldc);
