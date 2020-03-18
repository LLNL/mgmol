// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "entropy.h"
#include "MPIdata.h"
#include "fermi.h"
#include <cassert>
#include <cmath>
#include <iostream>

/* Evaluate the entropy function from occupations, f:
   s(f) = -f ln(f) - (1-f) ln(1-f).
   Returns the entropy = trace[s(f)] = sum {s(f)}_i
*/
template <typename T>
double entropy_eval(
    const std::vector<T>& f, std::vector<T>& s, const double occ_factor)
{
    const int dim = f.size();
    assert(s.size() == static_cast<unsigned int>(dim));

    const double tol = 1.e-15;

#ifndef NDEBUG
    const double tol_interval = 1.e-6;
#endif

    double entropy = 0.;

    for (int st = 0; st < dim; st++)
    {
        const double fi = (double)f[st];
        assert(fi >= 0. - tol_interval);
        assert(fi <= 1. + tol_interval);
        if (fi < tol)
        {
            s[st] = (T)(1. - fi) * log(1. - fi);
        }
        else if (fi > 1. - tol)
        {
            s[st] = (T)(fi * log(fi));
        }
        else
        {
            s[st] = (T)(fi * log(fi) + (1. - fi) * log(1. - fi));
        }
        entropy += s[st];
    }

    return (-occ_factor) * entropy; // in units of kbt
}

// Evaluate entropy given 'energies'
template <typename T>
double entropy_evalFromEnergies(const double mu, const int max_occ,
    const double kBT, const std::vector<T>& energies, std::vector<T>& s,
    const double occ_factor)
{
    std::vector<T> f((int)energies.size(), 0.);

    // calculate occupations based on fermi-dirac function
    fermi_distribution(mu, max_occ, kBT, energies, f);

    // call entropy function
    double ent = entropy_eval(f, s, occ_factor);

    return ent;
}

template double entropy_eval(const std::vector<double>& f,
    std::vector<double>& s, const double occ_factor);
template double entropy_evalFromEnergies(const double mu, const int max_occ,
    const double kBT, const std::vector<double>& energies,
    std::vector<double>& s, const double occ_factor);
