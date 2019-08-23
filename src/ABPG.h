// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
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

template <class T>
class ABPG : public OrbitalsStepper<T>
{
    Hamiltonian<T>* hamiltonian_;
    ProjectedMatricesInterface* proj_matrices_;
    Energy<T>* energy_;
    MGmol<T>* mgmol_strategy_;

    std::ostream& os_;

    Mixing<T>* wf_mix_;
    static bool pbset_;

    static Timer abpg_tm_;
    static Timer abpg_nl_update_tm_;
    static Timer comp_res_tm_;
    static Timer update_states_tm_;

    void update_states(T& orbitals, T& res, T& work_orbitals,
        const double precond_factor, const bool accelerate);

public:
    ABPG(Hamiltonian<T>* hamiltonian, ProjectedMatricesInterface* proj_matrices,
        MGmol<T>* mgmol_strategy, std::ostream& os)
        : hamiltonian_(hamiltonian),
          proj_matrices_(proj_matrices),
          mgmol_strategy_(mgmol_strategy),
          os_(os),
          wf_mix_(0)
    {
    }

    ~ABPG()
    {
        if (wf_mix_ != 0)
        {
            delete wf_mix_;
            wf_mix_ = 0;
        }
    }

    void setup(T&);

    int update(T& orbitals, Ions& ions, const double precond_factor,
        const bool orthof, T& work_orbitals, const bool accelerate,
        const bool print_res, const double atol);

    static void printTimers(std::ostream& os);
};

template <typename T>
Timer ABPG<T>::abpg_tm_("abpg_line_min");
template <typename T>
Timer ABPG<T>::abpg_nl_update_tm_("abpg_nl_update");
template <typename T>
Timer ABPG<T>::comp_res_tm_("abpg_comp_residuals_st");
template <typename T>
Timer ABPG<T>::update_states_tm_("abpg_update_states");

template <typename T>
bool ABPG<T>::pbset_ = false;

#endif
