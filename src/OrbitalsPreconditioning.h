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

template <class T>
class OrbitalsPreconditioning
{
private:
#ifdef HAVE_MAGMA
    using memory_space_type = MemorySpace::Device;
#else
    using memory_space_type = MemorySpace::Host;
#endif

    Preconditioning<MGPRECONDTYPE>* precond_;
    pb::GridFuncVector<MGPRECONDTYPE, memory_space_type>* gfv_work_;

    pb::GridFuncVector<MGPRECONDTYPE, memory_space_type>* gfv_work2_;

    // coefficient for preconditioning
    double gamma_;

    bool is_set_;

    // timers
    static Timer precond_tm_;

    Map2Masks* map2masks_;

public:
    OrbitalsPreconditioning()
    {
        is_set_    = false;
        precond_   = nullptr;
        gfv_work_  = nullptr;
        gfv_work2_ = nullptr;
    };

    ~OrbitalsPreconditioning();

    void setup(T& orbitals, const short mg_levels, const short lap_type,
        MasksSet*, const std::shared_ptr<LocalizationRegions>&);
    void precond_mg(T& orbitals);
    void setGamma(const pb::Lap<ORBDTYPE>& lapOper, const Potentials& pot,
        const short mg_levels, ProjectedMatricesInterface* proj_matrices);
    static void printTimers(std::ostream& os);
};

template <class T>
Timer OrbitalsPreconditioning<T>::precond_tm_(
    "OrbitalsPreconditioning::precond");

#endif
