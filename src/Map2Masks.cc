// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Map2Masks.h"

Map2Masks::Map2Masks(MasksSet* masks, const std::vector<int>& gids)
{
    assert(masks != nullptr);

    // associate each gid with a mask
    for (auto& gid : gids)
    {
        GridMask* pmask = masks->get_pmask(gid);
        assert(pmask != nullptr);
        gid2mask_.insert(std::pair<int, GridMask*>(gid, pmask));
    }
}

template <typename ScalarType>
void Map2Masks::apply(pb::GridFunc<ScalarType>& gu,
    const std::vector<std::vector<int>>& gids, const short level,
    const int color) const
{
    if (color == -1) return; // no mask applied
    if (gid2mask_.empty()) return; // no mask applied

    assert(gids.size() > 0);

    const short subdivx = gids.size();

    const pb::Grid& lgrid(gu.grid());
    const short nghosts = lgrid.ghost_pt();
    const int dim0      = lgrid.dim(0) / subdivx;
    const int incx      = lgrid.inc(0);
    const int lnumpt    = dim0 * incx;

    for (short iloc = 0; iloc < subdivx; iloc++)
    {
        int st = gids[iloc][color];

        if (st != -1)
        {
            std::map<int, GridMask*>::const_iterator it = gid2mask_.find(st);
            assert(it != gid2mask_.end());
            (it->second)->apply(gu, level, iloc);
        }
        else
        {
            int offset = (nghosts + dim0 * iloc) * incx;
            assert(offset + lnumpt < static_cast<int>(lgrid.sizeg()));
            ScalarType* pu = gu.uu() + offset;
            memset(pu, 0, lnumpt * sizeof(ScalarType));
        }
    }
}

template void Map2Masks::apply(pb::GridFunc<float>&,
    const std::vector<std::vector<int>>&, const short, const int) const;
template void Map2Masks::apply(pb::GridFunc<double>&,
    const std::vector<std::vector<int>>&, const short, const int) const;
