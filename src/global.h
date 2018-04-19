// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef __global_h__
#define __global_h__

// this file should be included in every MGmol file
// to enable global definitions, macros, ...

//#include "mgmol_memory.h"

#ifdef USE_MP
typedef float ORBDTYPE;
#else
typedef double ORBDTYPE;
#endif

/* lmasktype sets the data type for the mask coeffs */
typedef ORBDTYPE lmasktype;
//typedef float lmasktype;

typedef double RHODTYPE;
//typedef float RHODTYPE;

typedef double MATDTYPE;

typedef float MGPRECONDTYPE;

typedef double POTDTYPE;
//typedef float POTDTYPE;

typedef ORBDTYPE KBPROJDTYPE;

typedef float POISSONPRECONDTYPE;

#endif
