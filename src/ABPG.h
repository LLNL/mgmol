// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ABPG_H
#define MGMOL_ABPG_H

#include "Energy.h"
#include "Hamiltonian.h"
#include "MGmol.h"
#include "Mixing.h"
#include "OrbitalsStepper.h"
#include "Timer.h"

#include <iostream>

class ProjectedMatricesInterface;

template <class OrbitalsType>
class ABPG : public OrbitalsStepper<OrbitalsType>
{
    Hamiltonian<OrbitalsType>* hamiltonian_;
    Energy<OrbitalsType>* energy_;
    MGmol<OrbitalsType>* mgmol_strategy_;

    std::ostream& os_;

    Mixing<OrbitalsType>* wf_mix_;
    static bool pbset_;

    static Timer abpg_tm_;
    static Timer abpg_nl_update_tm_;
    static Timer comp_res_tm_;
    static Timer update_states_tm_;

    void update_states(OrbitalsType& orbitals, OrbitalsType& res,
        OrbitalsType& work_orbitals, const double precond_factor,
        const bool accelerate);

public:
    ABPG(Hamiltonian<OrbitalsType>* hamiltonian,
        MGmol<OrbitalsType>* mgmol_strategy, std::ostream& os)
        : hamiltonian_(hamiltonian),
          mgmol_strategy_(mgmol_strategy),
          os_(os),
          wf_mix_(nullptr)
    {
    }

    ~ABPG() override
    {
        if (wf_mix_ != nullptr)
        {
            delete wf_mix_;
            wf_mix_ = nullptr;
        }
    }

    void setup(OrbitalsType&) override;

    int updateWF(OrbitalsType& orbitals, Ions& ions,
        const double precond_factor, const bool orthof,
        OrbitalsType& work_orbitals, const bool accelerate,
        const bool print_res, const double atol) override;

    void restartMixing() override
    {
        if (wf_mix_ != nullptr) wf_mix_->restart();
    }

    static void printTimers(std::ostream& os);
};

template <typename OrbitalsType>
Timer ABPG<OrbitalsType>::abpg_tm_("abpg_line_min");
template <typename OrbitalsType>
Timer ABPG<OrbitalsType>::abpg_nl_update_tm_("abpg_nl_update");
template <typename OrbitalsType>
Timer ABPG<OrbitalsType>::comp_res_tm_("abpg_comp_residuals_st");
template <typename OrbitalsType>
Timer ABPG<OrbitalsType>::update_states_tm_("abpg_update_states");

template <typename OrbitalsType>
bool ABPG<OrbitalsType>::pbset_ = false;

#endif
