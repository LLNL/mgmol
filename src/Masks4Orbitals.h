// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "GridMask.h"
#include "Orbitals.h"

#include <map>
class MasksSet;
class LocalizationRegions;

class Masks4Orbitals
{
private:
    MasksSet* masks_;
    MasksSet* corrmasks_;

    std::map<int,GridMask*> gid_to_mask_;     // map gid->mask;
    std::map<int,GridMask*> gid_to_corrmask_; // map gid->corrmask;

    void associateGids2Masks(const LocalizationRegions& lrs);

public:
    
    Masks4Orbitals(MasksSet* masks, MasksSet* corrmasks, const LocalizationRegions& lrs);
    ~Masks4Orbitals();

    short checkOverlap(const int gid1, const int gid2, const short level);
    
    Vector3D center(const int gid)
    {
        assert( gid_to_mask_.find(gid)!=gid_to_mask_.end() );
        
        return gid_to_mask_[gid]->center();
    }
    
    GridMask& getMask(const int gid)
    {
        assert( masks_!=0 );
        assert( gid_to_mask_.find(gid)!=gid_to_mask_.end() );
        
        return *(gid_to_mask_[gid]);
    }
    GridMask& getCorrMask(const int gid)
    {
        assert( corrmasks_!=0 );
        assert( gid_to_corrmask_.find(gid)!=gid_to_corrmask_.end() );
        
        return *(gid_to_corrmask_[gid]);
    }
};
