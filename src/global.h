// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_GLOBAL_H
#define MGMOL_GLOBAL_H

// this file should be included in every MGmol file
// to enable global definitions, macros, ...

// #include "mgmol_memory.h"

#ifdef MGMOL_USE_MIXEDP
typedef float ORBDTYPE;
#else
typedef double ORBDTYPE;
#endif

/* lmasktype sets the data type for the mask coeffs */
typedef ORBDTYPE lmasktype;

typedef double RHODTYPE;

typedef double MATDTYPE;

typedef double POTDTYPE;

typedef double KBPROJDTYPE;

#endif
