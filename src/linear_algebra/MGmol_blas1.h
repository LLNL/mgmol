// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_MYBLAS1_H
#define MGMOL_MYBLAS1_H

#include "fc_mangle.h"

#include <cmath>
#include <string.h>

#define MY_VERSION 0
#define EPSILON 1.e-12

#ifdef __cplusplus
extern "C"
{
#endif

    void DAXPY(const int* const, const double* const, const double* const,
        const int* const, double*, const int* const);
    void SAXPY(const int* const, const float* const, const float* const,
        const int* const, float*, const int* const);
    void DSCAL(
        const int* const, const double* const, double*, const int* const);
    void SSCAL(const int* const, const float* const, float*, const int* const);
    void DCOPY(const int* const, const double* const, const int* const, double*,
        const int* const);
    void SCOPY(const int* const, const float* const, const int* const, float*,
        const int* const);
    double DDOT(const int* const, const double* const, const int* const,
        const double* const, const int* const);
    float SDOT(const int* const, const float* const, const int* const,
        const float* const, const int* const);
    double DNRM2(const int* const, const double* const, const int* const);
    float SNRM2(const int* const, const float* const, const int* const);
    int IDAMAX(const int* const, const double* const, const int* const);
    int ISAMAX(const int* const, const float* const, const int* const);
    void DROT(int*, double*, int*, double*, int*, double*, double*);
    void SROT(int*, float*, int*, float*, int*, float*, float*);
#ifdef __cplusplus
}
#endif

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

    DAXPY(&n, &alpha, a, &ione, b, &ione);
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

    DSCAL(&n, &alpha, a, &ione);
#endif
}

#endif
