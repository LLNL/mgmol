// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ReplicatedMatrix2SquareLocalMatrices_H
#define MGMOL_ReplicatedMatrix2SquareLocalMatrices_H

#include "ReplicatedMatrix.h"
#include "SquareLocalMatrices.h"
#include "Timer.h"

#include <memory>
#include <vector>

class ReplicatedMatrix2SquareLocalMatrices
{
    static ReplicatedMatrix2SquareLocalMatrices* pinstance_;

    static Timer convert_tm_;

    static std::vector<std::vector<int>> global_indexes_;

public:
    static ReplicatedMatrix2SquareLocalMatrices* instance()
    {
        if (pinstance_ == nullptr)
        {
            pinstance_ = new ReplicatedMatrix2SquareLocalMatrices();
        }
        return pinstance_;
    }

    ReplicatedMatrix2SquareLocalMatrices() { }

    static void setup(const std::vector<std::vector<int>>& gids)
    {
        global_indexes_ = gids;
    }

    ~ReplicatedMatrix2SquareLocalMatrices() { }

    void convert(const ReplicatedMatrix& dmat,
        SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& lmat);

    static void printTimers(std::ostream& os) { convert_tm_.print(os); }
};

#endif
