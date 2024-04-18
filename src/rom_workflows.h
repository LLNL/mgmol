// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef ROM_WORKFLOWS_H
#define ROM_WORKFLOWS_H

#include <cassert>
#include <iostream>
#include <iterator>
#include <vector>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_CNR
#include <mkl.h>
#endif

#include <mpi.h>

#include "Control.h"
#include "DistMatrix.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "MatricesBlacsContext.h"
#include "Mesh.h"
#include "PackedCommunicationBuffer.h"
#include "ReplicatedWorkSpace.h"
#include "SparseDistMatrix.h"
#include "magma_singleton.h"
#include "tools.h"

#include <fenv.h>
#include <sys/cdefs.h>
#include <time.h>

#include "OptionDescription.h"

#include "librom.h"

void readRestartFiles(po::variables_map &vm);

#endif  // ROM_WORKFLOWS_H
