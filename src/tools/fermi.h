// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef FERMI_H
#define FERMI_H

#include <vector>

double fermi_distribution(const double mu,
                          const int max_numst,
                          const double kBT,
                          const std::vector<double>& energies,
                          std::vector<double>& occ );

double fermi_distribution(const double mu,
                          const int max_numst,
                          const double kBT,
                          const std::vector<float>& energies,
                          std::vector<float>& occ );
#endif
