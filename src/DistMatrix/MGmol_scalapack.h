// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef MGMOL_MGMOL_USE_SCALAPACK_H
#define MGMOL_MGMOL_USE_SCALAPACK_H

#include "scalapack_mangle.h"

typedef const char* const Pchar;
typedef const int* const Pint;
typedef const double* const Pdouble;
typedef const float* const Pfloat;

extern "C"
{
#ifdef MGMOL_USE_SCALAPACK
    // PBLAS
    void pdswap(
        Pint, double*, Pint, Pint, Pint, Pint, double*, Pint, Pint, Pint, Pint);
    void psswap(
        Pint, float*, Pint, Pint, Pint, Pint, float*, Pint, Pint, Pint, Pint);
    void pdsymm(Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pint, Pint,
        Pdouble, Pint, Pint, Pint, Pdouble, double* const, Pint, Pint, Pint);
    void pssymm(Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint, Pint, Pint,
        Pfloat, Pint, Pint, Pint, Pfloat, float* const, Pint, Pint, Pint);
    void pdgemm(Pchar, Pchar, Pint, Pint, Pint, Pdouble, Pdouble, Pint, Pint,
        Pint, Pdouble, Pint, Pint, Pint, Pdouble, double* const, Pint, Pint,
        Pint);
    void pdgemv(Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pint, Pint, Pdouble,
        Pint, Pint, Pint, Pint, Pdouble, double* const, Pint, Pint, Pint, Pint);
    void pdsymv(Pchar, Pint, Pdouble, Pdouble, Pint, Pint, Pint, Pdouble, Pint,
        Pint, Pint, Pint, Pdouble, double* const, Pint, Pint, Pint, Pint);
    void psgemm(Pchar, Pchar, Pint, Pint, Pint, Pfloat, Pfloat, Pint, Pint,
        Pint, Pfloat, Pint, Pint, Pint, Pfloat, float* const, Pint, Pint, Pint);
    void psgemv(Pchar, Pint, Pint, Pfloat, Pfloat, Pint, Pint, Pint, Pfloat,
        Pint, Pint, Pint, Pint, Pfloat, float* const, Pint, Pint, Pint, Pint);
    void pssymv(Pchar, Pint, Pfloat, Pfloat, Pint, Pint, Pint, Pfloat, Pint,
        Pint, Pint, Pint, Pfloat, float* const, Pint, Pint, Pint, Pint);
    void pdsyrk(Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pint, Pint,
        Pdouble, double* const, Pint, Pint, Pint);
    void pssyrk(Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint, Pint, Pint,
        Pfloat, float* const, Pint, Pint, Pint);
    void pdtran(Pint, Pint, Pdouble, Pdouble, Pint, Pint, Pint, Pdouble,
        double* const, Pint, Pint, Pint);
    void pstran(Pint, Pint, Pfloat, Pfloat, Pint, Pint, Pint, Pfloat,
        float* const, Pint, Pint, Pint);
    void pdtrmm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint,
        Pint, Pint, double*, Pint, Pint, Pint);
    void pstrmm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint,
        Pint, Pint, float*, Pint, Pint, Pint);
    void pdtrsm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint,
        Pint, Pint, double*, Pint, Pint, Pint);
    void pstrsm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint,
        Pint, Pint, float*, Pint, Pint, Pint);
    void pdamax(Pint, double*, int*, double*, Pint, Pint, Pint, Pint);
    void psamax(Pint, float*, int*, float*, Pint, Pint, Pint, Pint);

    // MGMOL_USE_SCALAPACK
    void pdelset(double*, Pint, Pint, int*, Pdouble);
    void pselset(float*, Pint, Pint, int*, Pfloat);
    float pselget(Pchar, Pchar, float*, Pfloat, Pint, Pint, Pint);
    double pdelget(Pchar, Pchar, double*, Pdouble, Pint, Pint, Pint);
    double pdlatra(Pint, Pdouble, Pint, Pint, Pint);
    double pslatra(Pint, Pfloat, Pint, Pint, Pint);
    double pdlaset(
        char*, int*, int*, double*, double*, double*, int*, int*, int*);
    float pslaset(char*, int*, int*, float*, float*, float*, int*, int*, int*);
    // note: values of data pointed to by 3rd and 5th arguments of p?lacp3
    // may not actually be const, depending on the value of the last argument
    double pdlacp3(Pint, Pint, const double* const, Pint, const double* const,
        Pint, Pint, Pint, Pint);
    float pslacp3(Pint, Pint, const float* const, Pint, const float* const,
        Pint, Pint, Pint, Pint);
    double pdlange(Pchar, int*, int*, double*, int*, int*, int*, double*);
    float pslange(Pchar, int*, int*, float*, int*, int*, int*, float*);
    void pdtrtrs(Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pint, Pint, Pint,
        double*, Pint, Pint, Pint, int*);
    void pstrtrs(Pchar, Pchar, Pchar, Pint, Pint, Pfloat, Pint, Pint, Pint,
        float*, Pint, Pint, Pint, int*);
    void pdgemr2d(Pint, Pint, Pdouble, Pint, Pint, Pint, double* const, Pint,
        Pint, Pint, Pint);
    void psgemr2d(Pint, Pint, Pfloat, Pint, Pint, Pint, float* const, Pint,
        Pint, Pint, Pint);
    void pigemr2d(
        Pint, Pint, Pint, Pint, Pint, Pint, int*, Pint, Pint, Pint, Pint);
    void pdpotrf(Pchar, int*, double*, int*, int*, int*, int*);
    void pspotrf(Pchar, int*, float*, int*, int*, int*, int*);
    void pdpotrs(Pchar, int*, int*, double*, int*, int*, int*, double*, int*,
        int*, int*, int*);
    void pspotrs(Pchar, int*, int*, float*, int*, int*, int*, float*, int*,
        int*, int*, int*);
    void pdgetrf(int*, int*, double*, int*, int*, int*, int*, int*);
    void psgetrf(int*, int*, float*, int*, int*, int*, int*, int*);
    void pdgetrs(Pchar, int*, int*, double*, int*, int*, int*, int*, double*,
        int*, int*, int*, int*);
    void psgetrs(Pchar, int*, int*, float*, int*, int*, int*, int*, float*,
        int*, int*, int*, int*);
    void pdpotri(Pchar, int*, double*, int*, int*, int*, int*);
    void pspotri(Pchar, int*, float*, int*, int*, int*, int*);
    void pdtrtri(Pchar, Pchar, int*, double*, int*, int*, int*, int*);
    void pdpocon(Pchar, int*, double*, int*, int*, int*, double*, double*,
        double*, int*, int*, int*, int*);
    void pspocon(Pchar, int*, float*, int*, int*, int*, float*, float*, float*,
        int*, int*, int*, int*);
    void pdsygst(Pint, Pchar, Pint, double*, Pint, Pint, Pint, Pdouble, Pint,
        Pint, Pint, double*, int*);
    void pssygst(Pint, Pchar, Pint, float*, Pint, Pint, Pint, Pfloat, Pint,
        Pint, Pint, float*, int*);
    void pdsyev(Pchar, Pchar, int*, double*, int*, int*, int*, double*, double*,
        int*, int*, int*, double*, int*, int*);
    void pssyev(Pchar, Pchar, int*, float*, int*, int*, int*, float*, float*,
        int*, int*, int*, float*, int*, int*);
    void pdgesvd(Pchar, Pchar, int*, int*, double*, int*, int*, int*, double*,
        double*, int*, int*, int*, double*, int*, int*, int*, double*, int*,
        int*);
    void psgesvd(Pchar, Pchar, int*, int*, float*, int*, int*, int*, float*,
        float*, int*, int*, int*, float*, int*, int*, int*, float*, int*, int*);
    // MGMOL_USE_SCALAPACK TOOLS
    int NUMROC(Pint, Pint, Pint, Pint, Pint);
    int INDXL2G(Pint, Pint, Pint, Pint, Pint);
    int INDXG2L(Pint, Pint, Pint, Pint, Pint);
#endif
}

#endif
