// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef SinCosOps_H
#define SinCosOps_H

#include "Timer.h"
#include "VariableSizeMatrix.h"

class LocGridOrbitals;

#include <vector>

class SinCosOps
{
private :
    static Timer compute_tm_;

public:

    static void compute(const LocGridOrbitals& orbitals,
                        std::vector<std::vector<double> >& a);
    static void computeSquare(const LocGridOrbitals& orbitals,
                        std::vector<std::vector<double> >& a);
    static void compute1D(const LocGridOrbitals& orbitals,
                          std::vector<std::vector<double> >& a,
                          const int dim_index);
    static void computeSquare1D(const LocGridOrbitals& orbitals,
                                std::vector<std::vector<double> >& a,
                                const int dim_index);
    static void compute2states(const LocGridOrbitals& orbitals,
                               std::vector<std::vector<double> >& a,
                               const int st1, const int st2);
    static void computeDiag2states(const LocGridOrbitals& orbitals,
                                   std::vector<std::vector<double> >& a,
                                   const int st1, const int st2);
    static void compute(const LocGridOrbitals& orbitals1,
                        const LocGridOrbitals& orbitals2,
                        std::vector<std::vector<double> >& a);
    static void computeDiag(const LocGridOrbitals& orbitals,
                            VariableSizeMatrix<sparserow>& mat,
                            const bool normalized_functions);

    static void printTimers(std::ostream& os)
    {
        compute_tm_.print(os);
    }
};

#endif

