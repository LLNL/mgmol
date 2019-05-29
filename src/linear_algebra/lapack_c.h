// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_LAPACK_H
#define MGMOL_LAPACK_H

#include "fc_mangle.h"

typedef const char* const Pchar;

#define dsyev DSYEV
#define dpotri DPOTRI
#define dpotrf DPOTRF
#define dpotrs DPOTRS
#define dsygst DSYGST
#define dtrtrs DTRTRS
#define dpocon DPOCON
#define dsygv DSYGV
#define dlange DLANGE

extern "C"
{
    void dsyev(Pchar, Pchar, const int* const, double*, const int* const,
        double*, double*, const int* const, int*);
    void dsygv(const int* const, Pchar, Pchar, const int* const, double*,
        const int* const, double*, const int* const, double*, double*,
        const int* const, int*);
    void dpotri(Pchar, const int* const, double*, const int* const, int*);
    void dpotrf(Pchar, const int* const, double*, const int* const, int*);
    void dpotrs(Pchar, const int* const, const int* const, double*,
        const int* const, double*, const int* const, int*);
    void dgetrf(int*, int*, double*, int*, int*, int*);
    void dgetrs(Pchar, int*, int*, double*, int*, int*, double*, int*, int*);
    void dpocon(Pchar, const int* const, double*, const int* const, double*,
        double*, double*, const int* const, int*);
    void dtrtrs(Pchar, Pchar, Pchar, const int* const, const int* const,
        double*, const int* const, double*, const int* const, int*);
    void dtrtri(Pchar, Pchar, int*, double*, int*, int*);
    void dsygst(const int* const, Pchar, const int* const, double*,
        const int* const, double*, const int* const, int*);
    void dgesvd(Pchar, Pchar, int*, int*, double*, int*, double*, double*, int*,
        double*, int*, double*, int*, int*);
    double dlange(Pchar, int*, int*, double*, int*, double*);
}

#endif
