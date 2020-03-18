// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_SINCOSOPS_H
#define MGMOL_SINCOSOPS_H

#include "Timer.h"
#include "VariableSizeMatrix.h"

#include <vector>

template <class T>
class SinCosOps
{
private:
    static Timer compute_tm_;

public:
    static void compute(const T& orbitals, std::vector<std::vector<double>>& a);
    static void computeSquare(
        const T& orbitals, std::vector<std::vector<double>>& a);
    static void compute1D(const T& orbitals,
        std::vector<std::vector<double>>& a, const int dim_index);
    static void computeSquare1D(const T& orbitals,
        std::vector<std::vector<double>>& a, const int dim_index);
    static void compute2states(const T& orbitals,
        std::vector<std::vector<double>>& a, const int st1, const int st2);
    static void computeDiag2states(const T& orbitals,
        std::vector<std::vector<double>>& a, const int st1, const int st2);
    static void compute(const T& orbitals1, const T& orbitals2,
        std::vector<std::vector<double>>& a);
    static void computeDiag(const T& orbitals,
        VariableSizeMatrix<sparserow>& mat, const bool normalized_functions);

    static void printTimers(std::ostream& os) { compute_tm_.print(os); }
};

template <class T>
Timer SinCosOps<T>::compute_tm_("SinCosOps::compute_tm");

#endif
