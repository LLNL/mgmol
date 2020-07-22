// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// Class to allocate orbitals into "global" functions
#ifndef MGMOL_FUNCTIONSPACKING_H
#define MGMOL_FUNCTIONSPACKING_H

#include "SymmetricMatrix.h"

#include <list>
#include <map>
#include <memory>
#include <set>
#include <vector>

class LocalizationRegions;

class FunctionsPacking
{
private:
    std::map<int, short> gid2color_; // tells where one gid is allocated
    int global_size_;
    short num_colors_;

    MPI_Comm comm_;

    void setup(std::shared_ptr<LocalizationRegions> lrs, const bool global);

    void getColors(const SymmetricMatrix<char>& overlaps,
        std::list<std::list<int>>& colors);
    short checkOverlapLRs(std::shared_ptr<LocalizationRegions> lrs,
        const int gid1, const int gid2) const;
    void initOrbiOverlapLocal(std::shared_ptr<LocalizationRegions> lrs,
        const short level, SymmetricMatrix<char>& orbi_overlap);
    void initOrbiOverlapGlobal(std::shared_ptr<LocalizationRegions> lrs,
        const short level, SymmetricMatrix<char>& orbi_overlap);

public:
    FunctionsPacking(std::shared_ptr<LocalizationRegions> lrs,
        const bool global, const MPI_Comm comm);
    FunctionsPacking(const FunctionsPacking&);

    // return color of gid if exists locally, otherwise return -1
    short getColor(const int gid) const
    {
        assert(gid >= 0);

        std::map<int, short>::const_iterator it = gid2color_.find(gid);
        if (it == gid2color_.end())
            return -1;
        else
            return it->second;
    }

    int chromatic_number() const { return num_colors_; }
};

#endif
