// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// Adapted from MPIdata.h,v 1.5 2002/01/10 00:36:35 fgygi

#ifndef MPIDATA_H
#define MPIDATA_H

#ifdef USE_MPI

#include <mpi.h>
#include <iostream>
#include <fstream>
namespace MPIdata
{
  extern int mype;     // rank of this process
  extern bool onpe0;
  extern std::ostream* sout;
  extern std::ostream* serr;
};  
using namespace MPIdata;

#else

  const int mype = 0;
  const bool onpe0 = true;
  std::ostream* sout = &std::cout;
  std::ostream* serr = &std::cerr;

  typedef int MPI_Comm;

#endif

#endif
