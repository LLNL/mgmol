// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "SymmetricMatrix.h"

#include <iostream>
#include <list>

void colorRLF(const SymmetricMatrix<char>& overlaps,
    std::list<std::list<int>>& colored_gids, const bool, std::ostream&);

void greedyColor(const SymmetricMatrix<char>& overlaps,
    std::list<std::list<int>>& colored_gids, const bool, std::ostream&);

void greedyMC(int n, int* ja, int* ia, int* num_colors, int* iord, int* colors);
