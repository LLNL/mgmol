// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef LAPACK_H
#define LAPACK_H

typedef const char* const Pchar;

#ifdef UPCASE

  #define               dposv   SPOSV
  #define               dsyev   SSYEV
  #define               dgeev   SGEEV
  #define               dpotri  SPOTRI
  #define               dpotrf  SPOTRF
  #define               dpotrs  SPOTRS
  #define               dsygst  SSYGST
  #define               dtrtrs  STRTRS
  #define               dpocon  SPOCON
  #define               dsygv   SSYGV

#endif

#ifdef ADD_

  #define               dposv   dposv_
  #define               dsyev   dsyev_
  #define               dgeev   dgeev_
  #define               dpotri  dpotri_
  #define               dpotrf  dpotrf_
  #define               dpotrs  dpotrs_
  #define               dpocon  dpocon_
  #define               dsygst  dsygst_
  #define               dtrtrs  dtrtrs_
  #define               dpocon  dpocon_
  #define               dsygv   dsygv_

#endif



extern "C"{

void dposv(Pchar, const int* const, const int* const, double *, 
           const int* const, double *, const int* const, int *);
void dsyev(Pchar, Pchar, const int* const, double *, const int* const, 
           double *, double *, const int* const, int*);
void dgeev(Pchar, Pchar, const int* const, double*, const int* const, double*,
           double*, double*, const int* const, double*, const int* const, 
           double*, const int* const, int*);
void dsygv(const int* const, Pchar, Pchar, const int* const, double*,
           const int* const, double*,
           const int* const, double*, double*, const int* const, int*);
void dpotri(Pchar, const int* const, double*, const int* const, int*);
void dpotrf(Pchar, const int* const, double*, const int* const, int*);
void dpotrs(Pchar, const int* const, const int* const, double*, 
            const int* const, double*, const int* const, int*);
void dpocon(Pchar, const int* const, double *, const int* const, double *,
              double *, double *, const int* const,
             int *);
void dtrtrs(Pchar, Pchar, Pchar, const int* const, const int* const, double*, 
            const int* const, double*, const int* const, int*);
void dsygst(const int* const, Pchar, const int* const, double*, 
            const int* const, double*, const int* const, int*);


}

#endif
