// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_BLAS3_H
#define MGMOL_BLAS3_H

#include "fc_mangle.h"

typedef const char* const Pchar;
typedef const int* const Pint;
typedef const double* const Pdouble;
typedef const float* const Pfloat;

#define dgemm DGEMM
#define sgemm SGEMM
#define dsymm DSYMM
#define dgesv DGESV
#define dsyrk DSYRK
#define ssyrk SSYRK
#define dtrmm DTRMM
#define dtrsm DTRSM
#define strsm STRSM


extern "C"
{
    void dsyrk(Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pdouble,
        double* const, Pint);
    void ssyrk(Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint, Pfloat,
        float* const, Pint);
    void dgemm(Pchar, Pchar, Pint, Pint, Pint, Pdouble, Pdouble, Pint, Pdouble,
        Pint, Pdouble, double* const, Pint);
    void sgemm(Pchar, Pchar, Pint, Pint, Pint, Pfloat, Pfloat, Pint, Pfloat,
        Pint, Pfloat, float* const, Pint);
    void dsymm(Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pdouble, Pint,
        Pdouble, double* const, Pint);
    void dtrmm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint,
        double* const, Pint);
    void dtrsm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint,
        double* const, Pint);
    void strsm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint,
        float* const, Pint);
}

#endif
