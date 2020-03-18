// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_PRECONDITIONING_H
#define MGMOL_PRECONDITIONING_H

#include "GridFunc.h"
#include "GridMask.h"
#include "Lap.h"
#include "tools.h"

#include <map>
#include <vector>

template <typename T>
class Preconditioning
{
    short max_levels_;
    short lap_type_;
    short bc_[3];

    std::vector<pb::Lap<T>*> precond_oper_;
    std::vector<pb::Grid*> grid_;
    std::vector<pb::GridFunc<T>*> gf_work_;
    std::vector<pb::GridFunc<T>*> gf_rcoarse_;
    std::vector<pb::GridFunc<T>*> gf_newv_;

    std::vector<pb::GridFuncVector<T>*> gfv_work_;
    std::vector<pb::GridFuncVector<T>*> gfv_rcoarse_;
    std::vector<pb::GridFuncVector<T>*> gfv_newv_;

    std::map<int, GridMask*> gid2mask_; // map state->mask;
    std::vector<std::vector<int>> overlapping_gids_;

    void app_mask(
        pb::GridFunc<T>&, const short level, const int istate = -1) const;
    void app_mask(pb::GridFuncVector<T>&, const short level) const;

public:
    Preconditioning(const short lap_type, const short maxlevels,
        const pb::Grid& grid, const short bc[3]);

    Preconditioning(const Preconditioning&);

    ~Preconditioning();
    void clear();

    void setup(std::map<int, GridMask*>& st_to_mask,
        const std::vector<std::vector<int>>& gid);
    void reset(std::map<int, GridMask*>& st_to_mask,
        const std::vector<std::vector<int>>& gid);
    void mg(pb::GridFunc<T>& gf_v, const pb::GridFunc<T>& gf_f,
        const short level, const int istate = -1);
    void mg(pb::GridFuncVector<T>& gf_v, const pb::GridFuncVector<T>& gf_f,
        const short level);
    short max_levels() const { return max_levels_; }
};
#endif
