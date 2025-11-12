// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_OrbitalsPreconditioning_H
#define MGMOL_OrbitalsPreconditioning_H

#include "GridFuncVector.h"
#include "Lap.h"
#include "Map2Masks.h"
#include "Preconditioning.h"

#include <memory>

class Masks4Orbitals;
class MasksSet;
class ProjectedMatricesInterface;
class Potentials;
class LocalizationRegions;

template <class OrbitalsType, typename PDataType>
class OrbitalsPreconditioning
{
private:
#ifdef HAVE_MAGMA
    using memory_space_type = MemorySpace::Device;
#else
    using memory_space_type = MemorySpace::Host;
#endif

    std::shared_ptr<Preconditioning<PDataType>> precond_;

    // work arrays with preconditioner precision
    std::shared_ptr<pb::GridFuncVector<PDataType, memory_space_type>>
        gfv_work1_;
    std::shared_ptr<pb::GridFuncVector<PDataType, memory_space_type>>
        gfv_work2_;

    // tmp work array for case ORBDTYPE!=PDataType
    std::shared_ptr<pb::GridFuncVector<ORBDTYPE, memory_space_type>> gfv_work3_;

    short lap_type_;

    // coefficient for preconditioning
    double gamma_;

    bool is_set_;

    // timers
    static Timer precond_tm_;

    std::shared_ptr<Map2Masks> map2masks_;

public:
    OrbitalsPreconditioning() { is_set_ = false; };

    ~OrbitalsPreconditioning();

    void setup(OrbitalsType& orbitals, const short mg_levels,
        const short lap_type, MasksSet*,
        const std::shared_ptr<LocalizationRegions>&);
    void precond_mg(OrbitalsType& orbitals);
    void setGamma(const pb::Lap<ORBDTYPE>& lapOper, const Potentials& pot,
        const short mg_levels, ProjectedMatricesInterface* proj_matrices);
    static void printTimers(std::ostream& os);
};

template <class OrbitalsType, typename PDataType>
Timer OrbitalsPreconditioning<OrbitalsType, PDataType>::precond_tm_(
    "OrbitalsPreconditioning::precond");

#endif
