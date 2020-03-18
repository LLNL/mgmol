// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_BLAS2_H
#define MGMOL_BLAS2_H

typedef const char* const Pchar;
typedef const int* const Pint;
typedef const double* const Pdouble;
typedef const float* const Pfloat;

extern "C"
{
    void DSYMV(Pchar, Pint, Pdouble, Pdouble, Pint, Pdouble, Pint, Pdouble,
        double*, Pint);
    void DGEMV(Pchar, Pint, Pint, Pdouble, Pdouble, Pint, Pdouble, Pint,
        Pdouble, double*, Pint);
    void SGEMV(Pchar, Pint, Pint, Pfloat, Pfloat, Pint, Pfloat, Pint, Pfloat,
        float*, Pint);
}

#endif
