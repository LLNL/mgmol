// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Masks4Orbitals.h"
#include "Control.h"
#include "MasksSet.h"
using namespace std;

Masks4Orbitals::Masks4Orbitals(
    MasksSet* masks, MasksSet* corrmasks, const vector<int>& overlap_gids)
    : masks_(masks), corrmasks_(corrmasks)
{
    associateGids2Masks(overlap_gids);
}

Masks4Orbitals::~Masks4Orbitals() { }

void Masks4Orbitals::associateGids2Masks(const vector<int>& overlap_gids)
{
    Control& ct      = *(Control::instance());
    const int nmasks = (int)(masks_->size());
    if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout) << "Masks4Orbitals::associateGids2Masks()" << endl;
    }

    gid_to_mask_.clear();
    gid_to_corrmask_.clear();

    for (vector<int>::const_iterator it = overlap_gids.begin();
         it != overlap_gids.end(); it++)
    {
        int gid         = (*it);
        GridMask* maski = masks_->get_pmask(gid);
        assert(maski != nullptr);
        gid_to_mask_.insert(pair<int, GridMask*>(gid, maski));
    }

    if (corrmasks_ != nullptr)
        for (vector<int>::const_iterator it = overlap_gids.begin();
             it != overlap_gids.end(); it++)
        {
            int gid          = (*it);
            GridMask* maskci = corrmasks_->get_pmask(gid);
            assert(maskci != nullptr);
            gid_to_corrmask_.insert(pair<int, GridMask*>(gid, maskci));
        }

    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << "Initialize LocGridOrbitals with " << nmasks
                         << " masks" << endl;
}

short Masks4Orbitals::checkOverlap(
    const int gid1, const int gid2, const short level)
{
    map<int, GridMask*>::const_iterator it1 = gid_to_mask_.find(gid1);
    if (it1 == gid_to_mask_.end()) return 0;

    map<int, GridMask*>::const_iterator it2 = gid_to_mask_.find(gid2);
    if (it2 == gid_to_mask_.end()) return 0;

    GridMask* gm1 = it1->second;
    GridMask* gm2 = it2->second;
    assert(gm1 != nullptr);
    assert(gm2 != nullptr);

    return gm1->overlap_on_pe(*gm2, level);
}
