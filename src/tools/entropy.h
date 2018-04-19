// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef ENTROPY_H
#define ENTROPY_H

#include <vector>

double entropy_eval(const std::vector<double>& f, std::vector<double>& s, const double occ_factor);

double entropy_eval(const std::vector<float>& f, std::vector<float>& s, const double occ_factor);

double entropy_evalFromEnergies(const double mu,
                          const int max_occ,
                          const double kBT,
                          const std::vector<double>& energies, 
                          std::vector<double>& s, 
                          const double occ_factor);

double entropy_evalFromEnergies(const double mu,
                          const int max_occ,
                          const double kBT,
                          const std::vector<float>& energies, 
                          std::vector<float>& s, 
                          const double occ_factor);                       
                          
#endif
