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
#include "GridFuncVector.h"
#include "GridMask.h"
#include "Lap.h"
#include "tools.h"

#include <map>
#include <vector>

template <typename T>
class Preconditioning
{
public:
#ifdef HAVE_MAGMA
    using memory_space_type = MemorySpace::Device;
#else
    using memory_space_type = MemorySpace::Host;
#endif

private:
    short max_levels_;
    short lap_type_;
    short bc_[3];

    std::vector<pb::Lap<T>*> precond_oper_;
    std::vector<pb::Grid*> grid_;
    std::vector<pb::GridFunc<T>*> gf_work_;
    std::vector<pb::GridFunc<T>*> gf_rcoarse_;
    std::vector<pb::GridFunc<T>*> gf_newv_;

    std::vector<pb::GridFuncVector<T, memory_space_type>*> gfv_work_;
    std::vector<pb::GridFuncVector<T, memory_space_type>*> gfv_rcoarse_;
    std::vector<pb::GridFuncVector<T, memory_space_type>*> gfv_newv_;

public:
    Preconditioning(const short lap_type, const short maxlevels,
        const pb::Grid& grid, const short bc[3]);

    Preconditioning(const Preconditioning&);

    ~Preconditioning();
    void clear();

    void setup(const std::vector<std::vector<int>>& overlapping_gids);
    void reset(const std::vector<std::vector<int>>& overlapping_gids);
    void mg(pb::GridFuncVector<T, memory_space_type>& gf_v,
        const pb::GridFuncVector<T, memory_space_type>& gf_f,
        const short level);
    short max_levels() const { return max_levels_; }
};
#endif
