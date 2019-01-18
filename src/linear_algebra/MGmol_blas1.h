// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: MGmol_blas1.h,v 1.12 2009/04/22 22:53:32 jeanluc Exp $
#ifndef PB_MYBLAS1_H
#define PB_MYBLAS1_H

#include <cmath>
#include <string.h>

#ifdef USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#define MY_VERSION 0
#define EPSILON 1.e-12

#ifdef ADD_ // LINUX, SGI, SUN

#define daxpy daxpy_
#define saxpy saxpy_
#define dscal dscal_
#define sscal sscal_
#define dcopy dcopy_
#define scopy scopy_
#define ddot ddot_
#define sdot sdot_
#define dnrm2 dnrm2_
#define snrm2 snrm2_
#define dswap dswap_
#define sswap sswap_
#define idamax idamax_
#define isamax isamax_

#endif

#ifdef UPCASE // CRAY_T3E

#define daxpy SAXPY
#define dscal SSCAL
#define dcopy SCOPY
#define ddot SDOT
#define dnrm2 SNRM2
#define dswap SSWAP
#define idamax IDAMAX
#define isamax ISAMAX

#endif

#ifdef __cplusplus
extern "C"
{
#endif

    void daxpy(const int* const, const double* const, const double* const,
        const int* const, double*, const int* const);
    void saxpy(const int* const, const float* const, const float* const,
        const int* const, float*, const int* const);
    void dscal(
        const int* const, const double* const, double*, const int* const);
    void sscal(const int* const, const float* const, float*, const int* const);
    void dcopy(const int* const, const double* const, const int* const, double*,
        const int* const);
    void scopy(const int* const, const float* const, const int* const, float*,
        const int* const);
    double ddot(const int* const, const double* const, const int* const,
        const double* const, const int* const);
    float sdot(const int* const, const float* const, const int* const,
        const float* const, const int* const);
    double dnrm2(const int* const, const double* const, const int* const);
    float snrm2(const int* const, const float* const, const int* const);
    void dswap(
        const int* const, double*, const int* const, double*, const int* const);
    void sswap(
        const int* const, float*, const int* const, float*, const int* const);
    int idamax(const int* const, const double* const, const int* const);
    int isamax(const int* const, const float* const, const int* const);
#ifdef __cplusplus
}
#endif

double my_ddot(int, const double* const, const double* const);
double my_pddot(int, const double* const, const double* const, MPI_Comm);

inline void my_dcopy(const int n, const double* const a, double* b)
{
    int ione = 1;
    dcopy(&n, a, &ione, b, &ione);
}

inline void my_daxpy(
    const int n, const double alpha, const double* const a, double* b)
{
#if MY_VERSION
    register int i;

    if (fabs(alpha - 1.) < EPSILON)
    {
        for (i = 0; i < n; i++)
            b[i] += a[i];
    }
    else
    {
        for (i = 0; i < n; i++)
            b[i] += alpha * a[i];
    }
#else
    int ione = 1;

    daxpy(&n, &alpha, a, &ione, b, &ione);
#endif
}

inline void my_dscal(const int n, const double alpha, double* a)
{
#if MY_VERSION
    register int i;

    for (i = 0; i < n; i++)
        a[i] *= alpha;
#else
    int ione = 1;

    dscal(&n, &alpha, a, &ione);
#endif
}

inline double my_dnrm2(int n, double* a)
{
    double dot = my_ddot(n, a, a);
    return sqrt(dot);
}

inline double my_pdnrm2(int n, double* a, MPI_Comm comm)
{
    double dot = my_pddot(n, a, a, comm);
    return sqrt(dot);
}

#endif
