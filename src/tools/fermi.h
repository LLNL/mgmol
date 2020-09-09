// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_FERMI_H
#define MGMOL_FERMI_H

#include <vector>

template <typename T>
double fermi_distribution(const double mu, const int max_numst,
    const double kBT, const std::vector<T>& energies, std::vector<T>& occ);

template <typename T>
double compute_chemical_potential_and_occupations(
    const std::vector<T>& energies, const double width, const double nel,
    const int max_numst, const bool onpe0, std::vector<T>& occ);

#endif
