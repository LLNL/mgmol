// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_Map2Masks
#define MGMOL_Map2Masks

#include "GridFunc.h"
#include "GridMask.h"
#include "MasksSet.h"

#include <map>
#include <vector>

class Map2Masks
{
private:
    std::map<int, GridMask*> gid2mask_;

public:
    Map2Masks(MasksSet* masks, const std::vector<int>& overlapping_gids);

    template <typename ScalarType>
    void apply(pb::GridFunc<ScalarType>& gu,
        const std::vector<std::vector<int>>&, const short level,
        const int istate) const;
};

#endif
