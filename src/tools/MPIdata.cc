// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

////////////////////////////////////////////////////////////////////////////////
//
// MPIdata.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id$
// Adapted from MPIdata.C,v 1.4 2002/01/10 00:36:34 fgygi

#include "MPIdata.h"

#ifdef USE_MPI
int MPIdata::mype;     // rank of this process
bool MPIdata::onpe0;
std::ostream* MPIdata::sout=&std::cout; //default value
std::ostream* MPIdata::serr=&std::cerr; //default value
#endif