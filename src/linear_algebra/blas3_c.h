// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
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

extern "C"
{
    void DSYRK(Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pdouble,
        double* const, Pint);
    void SSYRK(Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint, Pfloat,
        float* const, Pint);
    void DGEMM(Pchar, Pchar, Pint, Pint, Pint, Pdouble, Pdouble, Pint, Pdouble,
        Pint, Pdouble, double* const, Pint);
    void SGEMM(Pchar, Pchar, Pint, Pint, Pint, Pfloat, Pfloat, Pint, Pfloat,
        Pint, Pfloat, float* const, Pint);
    void DSYMM(Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pdouble, Pint,
        Pdouble, double* const, Pint);
    void DTRMM(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint,
        double* const, Pint);
    void DTRSM(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pdouble, Pdouble, Pint,
        double* const, Pint);
    void STRSM(Pchar, Pchar, Pchar, Pchar, Pint, Pint, Pfloat, Pfloat, Pint,
        float* const, Pint);
}

#endif
