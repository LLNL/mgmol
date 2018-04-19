// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef BLAS3_H
#define BLAS3_H

typedef const char* const  Pchar;
typedef const int* const     Pint;
typedef const double* const  Pdouble;
typedef const float* const  Pfloat;


#ifdef ADD_

  #define    dgemm    dgemm_
  #define    sgemm    sgemm_  
  #define    dsymm    dsymm_
  #define    dsyrk    dsyrk_
  #define    ssyrk    ssyrk_
  #define    dtrmm    dtrmm_
  #define    dtrsm    dtrsm_
  #define    strsm    strsm_  

#endif

#ifdef UPCASE

  #define    dgemm    SGEMM
  #define    dsymm    SSYMM
  #define    dgesv    SGESV
  #define    dsyrk    SSYRK
  #define    dtrmm    STRMM
  #define    dtrsm    STRSM

#endif


extern "C"{

void dsyrk(Pchar, Pchar, Pint, Pint,
           Pdouble, Pdouble, Pint, Pdouble, double* const, Pint);
void ssyrk(Pchar, Pchar, Pint, Pint,
           Pfloat, Pfloat, Pint, Pfloat, float* const, Pint);
void dgemm(Pchar, Pchar, 
           Pint, Pint, Pint,
           Pdouble, Pdouble, Pint, 
           Pdouble , Pint, Pdouble, double* const, Pint);
void sgemm(Pchar, Pchar, 
           Pint, Pint, Pint,
           Pfloat, Pfloat, Pint, 
           Pfloat , Pint, Pfloat, float* const, Pint);           
void dsymm(Pchar, Pchar, Pint, Pint,
           Pdouble, Pdouble, Pint, 
           Pdouble, Pint, Pdouble, double* const, Pint);
void dtrmm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, 
           Pdouble, Pdouble, Pint, double* const, Pint);
void dtrsm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, 
           Pdouble, Pdouble, Pint , double* const, Pint);
void strsm(Pchar, Pchar, Pchar, Pchar, Pint, Pint, 
           Pfloat, Pfloat, Pint , float* const, Pint);
           
}

#endif
