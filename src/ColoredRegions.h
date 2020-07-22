// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
//
#ifndef COLORED_REGIONS_H
#define COLORED_REGIONS_H

#include "FunctionsPacking.h"
#include "Vector3D.h"

#include <map>
#include <set>
#include <vector>

class LocalizationRegions;

typedef struct Sphere
{
    Vector3D center;
    double radii;
    int gid;
} Sphere;

class ColoredRegions
{
private:
    // orbitals centers (x,y,z) and radii for each global array
    // 1st index=color, 2nd index=center+radii+gid
    std::multimap<short, Sphere> loc_regions_;

    FunctionsPacking& pack_;

    void setup(std::shared_ptr<LocalizationRegions> lrs, const bool global);

public:
    ColoredRegions(FunctionsPacking& pack,
        std::shared_ptr<LocalizationRegions> lrs, const bool global);

    void getPossibleColors(const Vector3D& center, std::set<int>& colors)
    {
        colors.clear();
        std::multimap<short, Sphere>::iterator p = loc_regions_.begin();
        while (p != loc_regions_.end())
        {
            if (p->second.center == center) colors.insert(p->first);
            p++;
        }
    }

    // get localization centers stored in function color
    int getLocCentersAndRadii4color(
        const short color, std::vector<double>& centers_and_radii) const;

    int getAllCentersAndRadii4color(
        const short color, std::vector<double>& centers_and_radii) const;

    int getLocGids4color(const short color, std::vector<int>& data) const;

    int getAllGids4color(const short color, std::vector<int>& gids) const;
};

#endif
