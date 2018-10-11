// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
//
#ifndef COLORED_REGIONS_H
#define COLORED_REGIONS_H

#include "Vector3D.h"
#include "FunctionsPacking.h"

#include <map>
#include <vector>
#include <set>

class LocalizationRegions;

typedef struct Sphere
{
  Vector3D center;
  double   radii;
  int      gid;
} Sphere;


class ColoredRegions
{
private:
    // orbitals centers (x,y,z) and radii for each global array
    // 1st index=color, 2nd index=center+radii+gid
    std::multimap<short,Sphere>  loc_regions_;

    FunctionsPacking& pack_;

    void setup(LocalizationRegions& lrs,
               const bool global);

public:
    ColoredRegions(FunctionsPacking& pack,
                   LocalizationRegions& lrs, const bool global);

    void getPossibleColors(const Vector3D center, std::set<int>& colors)
    {
        colors.clear();
        std::multimap<short,Sphere>::iterator p=loc_regions_.begin();
        while( p!=loc_regions_.end() )
        {
            if( p->second.center==center )colors.insert(p->first);
            p++;
        }
    }

    // get localization centers stored in function color
    int getLocCentersAndRadii4color(const short color,
        std::vector<double>& centers_and_radii)const;

    int getAllCentersAndRadii4color(
        const short color, std::vector<double>& centers_and_radii)const;

    int getLocGids4color(const short color, std::vector<int>& data)const;

    int getAllGids4color(const short color, std::vector<int>& gids)const;

};

#endif

