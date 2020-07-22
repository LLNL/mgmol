// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_MASKSSET_H
#define MGMOL_MASKSSET_H

#include "GridMask.h"
#include "MPIdata.h"
#include "Vector3D.h"
#include "global.h"

#include <map>
#include <memory>

class LocalizationRegions;

class MasksSet
{

private:
    // containers for actual mask data
    std::map<int, GridMask*> pgrid_masks_;

    // type
    const bool type_corr_mask_;

    const short mg_levels_;

    MasksSet(const MasksSet&);

    void allocate(const std::vector<int>& gids);

public:
    MasksSet(const bool type_corr_mask, const short levels = 0)
        : type_corr_mask_(type_corr_mask), mg_levels_(levels)
    {
    }

    MasksSet(const std::shared_ptr<LocalizationRegions> lrs,
        const bool type_corr_mask, const short levels = 0,
        const double override_radius = 0.)
        : type_corr_mask_(type_corr_mask), mg_levels_(levels)
    {
        setup(lrs, override_radius);
    }
    MasksSet& operator=(const MasksSet& masks);

    // Destructor
    ~MasksSet() { clear(); }

    void clear();
    void setup(const std::shared_ptr<LocalizationRegions>,
        const double override_radius = 0.);
    void update(const std::shared_ptr<LocalizationRegions>);
    int initialize(const std::shared_ptr<LocalizationRegions>,
        const unsigned short, const double override_radius = 0.);
    GridMask* get_pmask(const int i)
    {
        assert(pgrid_masks_.size() > 0);
        std::map<int, GridMask*>::const_iterator it = pgrid_masks_.find(i);
        if (it == pgrid_masks_.end())
            return nullptr;
        else
            return it->second;
    }
    size_t size() const { return pgrid_masks_.size(); }
};

#endif
