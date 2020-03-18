// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ENTROPY_H
#define MGMOL_ENTROPY_H

#include <vector>

template <typename T>
double entropy_eval(
    const std::vector<T>& f, std::vector<T>& s, const double occ_factor);

template <typename T>
double entropy_evalFromEnergies(const double mu, const int max_occ,
    const double kBT, const std::vector<T>& energies, std::vector<T>& s,
    const double occ_factor);

#endif
