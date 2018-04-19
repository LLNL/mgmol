// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef BLAS2_H
#define BLAS2_H

typedef const char* const  Pchar;
typedef const int* const     Pint;
typedef const double* const  Pdouble;


#ifdef ADD_

#define    dsymv    dsymv_
#define    dgemv    dgemv_

#endif

#ifdef UPCASE

#define    dsymv    SSYMV
#define    dgemv    SGEMV

#endif


extern "C"{

void dsymv(Pchar, Pint, Pdouble, Pdouble, 
           Pint, Pdouble, Pint, Pdouble, double *, Pint);
void dgemv(Pchar, Pint, Pint, Pdouble, 
           Pdouble, Pint, Pdouble, Pint, Pdouble, double *, Pint);
}

#endif
