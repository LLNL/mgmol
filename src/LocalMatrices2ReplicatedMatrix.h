// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_LocalMatrices2ReplicatedMatrix_H
#define MGMOL_LocalMatrices2ReplicatedMatrix_H

#include "LocalMatrices.h"
#include "ReplicatedMatrix.h"
#include "Timer.h"

#include <vector>

// Add matrix elements corresponding to subdomains at their right place
// into a ReplicatedMatrix
// Important Note: Neglect contributions smaller than tol!
// (may lead to results dependent on number of CPUs)

class LocalMatrices2ReplicatedMatrix
{
private:
    static LocalMatrices2ReplicatedMatrix* pinstance_;

    static Timer convert_tm_;

    static std::vector<std::vector<int>> global_indexes_;

    static double tol_mat_elements;

public:
    static LocalMatrices2ReplicatedMatrix* instance()
    {
        if (pinstance_ == nullptr)
        {
            pinstance_ = new LocalMatrices2ReplicatedMatrix();
        }
        return pinstance_;
    }

    LocalMatrices2ReplicatedMatrix() { }

    static void setup(const std::vector<std::vector<int>>& gids)
    {
        global_indexes_ = gids;
    }

    void convert(const LocalMatrices<double, MemorySpace::Host>& src,
        ReplicatedMatrix& dst, const int numst,
        const double tol = tol_mat_elements) const;

    void accumulate(const LocalMatrices<double, MemorySpace::Host>& src,
        ReplicatedMatrix& dst, const double tol = tol_mat_elements) const;

    static void printTimers(std::ostream& os) { convert_tm_.print(os); }
};

#endif
