// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_PROTOTYPES_H
#define MGMOL_PROTOTYPES_H

#include "global.h"

class Ions;
class KBPsiMatrixSparse;

void get_vnlpsi(const Ions& ions, const std::vector<std::vector<int>>&,
    const int, const KBPsiMatrixSparse* const kbpsi, ORBDTYPE* const);
double getLAeigen(const double tol, const int maxit, Ions& ions);

#endif
