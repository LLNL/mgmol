// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// Adapted from MPIdata.h,v 1.5 2002/01/10 00:36:35 fgygi

#ifndef MPIDATA_H
#define MPIDATA_H

#include <fstream>
#include <iostream>
#include <mpi.h>
namespace MPIdata
{
extern int mype; // rank of this process
extern bool onpe0;
extern std::ostream* sout;
extern std::ostream* serr;
}
using namespace MPIdata;

#endif
