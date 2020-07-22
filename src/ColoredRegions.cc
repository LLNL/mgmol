// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
//

#include "ColoredRegions.h"
#include "Control.h"
#include "LocalizationRegions.h"

using namespace std;

ColoredRegions::ColoredRegions(FunctionsPacking& pack,
    std::shared_ptr<LocalizationRegions> lrs, const bool global)
    : pack_(pack)
{
    setup(lrs, global);
}

void ColoredRegions::setup(
    std::shared_ptr<LocalizationRegions> lrs, const bool global)
{
    std::vector<int> gids;
    if (global)
        lrs->getGidsGlobal(gids);
    else
        gids = lrs->getOverlapGids();

    loc_regions_.clear();
    Control& ct = *(Control::instance());

    if (ct.verbose > 0)
        printWithTimeStamp(" ColoredRegions::setLocRegions", cout);

    vector<int>::const_iterator ic = gids.begin();
    if (lrs->globalNumLRs() > 1)
    {
        // loop over gids
        while (ic != gids.end())
        {

            const int gid     = *ic;
            const short color = pack_.getColor(gid);

            assert(gid >= 0);
            assert(gid < (int)lrs->globalNumLRs());
            if (lrs->radius(gid) > 0.1) // radius=-1 if gid not know locally
            {
                Sphere lr;
                lr.center = lrs->getCenter(gid);
                lr.radii  = lrs->radius(gid);
                assert(lr.radii > 0.);
                lr.gid = gid;
                loc_regions_.insert(pair<short, Sphere>(color, lr));
            }
            ic++;
        }
    }
    else
    {
        assert(lrs->globalNumLRs() == 0 || lrs->globalNumLRs() == 1);

        short color = 0;
        while (ic != gids.end())
        {

            Sphere lr;
            lr.center = lrs->getCenter(0);
            lr.radii  = lrs->radius(0);
            lr.gid    = *ic;
            loc_regions_.insert(pair<short, Sphere>(color, lr));

            color++;
            ic++;
        }
    }
}

int ColoredRegions::getLocCentersAndRadii4color(
    const short color, vector<double>& data) const
{
    data.clear();
    multimap<short, Sphere>::const_iterator p = loc_regions_.find(color);
    if (p == loc_regions_.end()) return 0;

    while (p != loc_regions_.upper_bound(color))
    {
        double arr[] = { p->second.center[0], p->second.center[1],
            p->second.center[2], p->second.radii };
        data.insert(data.end(), arr, arr + 4);
        p++;
    }
    return data.size() / 4;
}

int ColoredRegions::getAllCentersAndRadii4color(
    const short color, vector<double>& centers_and_radii) const
{
    vector<double> local_centers_and_radii;
    getLocCentersAndRadii4color(color, local_centers_and_radii);

    centers_and_radii.clear();
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allGatherV(local_centers_and_radii, centers_and_radii);

    return centers_and_radii.size() / 4;
}

// get localization centers and radii stored in function color
int ColoredRegions::getLocGids4color(const short color, vector<int>& data) const
{
    data.clear();
    multimap<short, Sphere>::const_iterator p = loc_regions_.find(color);
    if (p == loc_regions_.end()) return 0;

    while (p != loc_regions_.upper_bound(color))
    {
        data.push_back(p->second.gid);
        p++;
    }
    return data.size();
}

int ColoredRegions::getAllGids4color(const short color, vector<int>& gids) const
{
    vector<int> lgids;
    getLocGids4color(color, lgids);

    gids.clear();
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allGatherV(lgids, gids);

    return gids.size();
}
