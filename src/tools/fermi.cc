// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "fermi.h"
#include <cassert>
#include <cmath>

using namespace std;

// fill occupations according to the fermi-dirac distribution:
//   f(x) = 1 / (1 + Exp[kB*(E-mu)])
// states 0 to max_occ-1, set to 0 states above max_occ-1
// returns the number of electrons
template <typename T>
double fermi_distribution(const double mu, const int max_occ, const double kBT,
    const vector<T>& energies, vector<T>& occ)
{
    assert(kBT >= 0.);
    assert(occ.size() == energies.size());

    const int dim = occ.size();

    for (int i = max_occ; i < dim; i++)
        occ[i] = 0.;

    // if(onpe0)(*MPIdata::sout)<<"fermi_distribution() with mu="
    //             <<mu<<" and kBT "<<kBT<<endl;
    double sum_occ = 0.;

    if (kBT > 1.e-10)
    {
        const double beta = 1. / kBT;
        for (int st = 0; st < max_occ; st++)
        {
            const double t1 = (energies[st] - mu) * beta;

            if (t1 > 0.)
            {
                const double t2 = exp(-t1);
                occ[st]         = t2 / (1. + t2);
            }
            else
            {
                const double t2 = exp(t1);
                occ[st]         = (T)(1. / (1. + t2));
            }

            sum_occ += occ[st];
        }
    }
    else // zero temperature
    {
        for (int st = 0; st < max_occ; st++)
        {
            const double t1 = (energies[st] - mu);

            if (t1 > 1.e-6)
            {
                occ[st] = 0.;
            }
            else
            {
                occ[st] = 1.;
            }
            sum_occ += occ[st];
        }
    }
    return sum_occ;
}

template double fermi_distribution(const double mu, const int max_occ,
    const double kBT, const vector<double>& energies, vector<double>& occ);
