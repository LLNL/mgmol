// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef MASKSSET_H
#define MASKSSET_H

#include "Vector3D.h"
#include "MPIdata.h"
#include "GridMask.h"
#include "global.h"

#include <map>

class LocalizationRegions;

class MasksSet{

private:

    // containers for actual mask data
    std::map<int,GridMask*> pgrid_masks_;

    // type
    const bool type_corr_mask_;

    const short   mg_levels_;

    MasksSet(const MasksSet&);
    
    void allocate(const std::vector<int>& gids);

public:
    
    MasksSet(const bool type_corr_mask,
             const short levels=0):
        type_corr_mask_(type_corr_mask),
        mg_levels_(levels)
    {
    }
    
    MasksSet(const LocalizationRegions& lrs,
             const bool type_corr_mask,
             const short levels=0,
             const double override_radius=0.):
        type_corr_mask_(type_corr_mask),
        mg_levels_(levels)
    {
        setup(lrs, override_radius);
    }
    MasksSet& operator=(const MasksSet& masks);
        
    // Destructor
    ~MasksSet()
    {
        clear();
    }
    
    void clear();
    void setup(const LocalizationRegions&,
               const double override_radius=0.);
    void update(const LocalizationRegions&);
    int initialize(const LocalizationRegions&, const unsigned short,
                   const double override_radius=0.);
    GridMask* get_pmask(const int i)
    {
        assert( pgrid_masks_.size()>0 );
        std::map<int,GridMask*>::const_iterator it=pgrid_masks_.find(i);
        if( it==pgrid_masks_.end() )return 0;
        else return it->second;
    }
    size_t size()const{ return pgrid_masks_.size(); }
};

#endif