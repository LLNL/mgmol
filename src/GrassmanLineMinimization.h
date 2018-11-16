// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef _GRASSMANLINEMINIMIZATION_H_
#define _GRASSMANLINEMINIMIZATION_H_

#include "Ions.h"
#include "LocGridOrbitals.h"
#include "OrbitalsStepper.h"
#include "ProjectedMatricesInterface.h"

#include <iostream>

class MGmol;

class GrassmanLineMinimization : public OrbitalsStepper
{
protected:
    Hamiltonian* hamiltonian_;
    ProjectedMatricesInterface* proj_matrices_;
    MGmol* mgmol_strategy_;
    std::ostream& os_;

    static bool accelerate_;
    static bool pbset_;
    static bool conjugate_;

    Ions* ptr2ions_;

    // Pointers to residual and search direction
    LocGridOrbitals* new_grad_;
    LocGridOrbitals* new_pcgrad_;

    LocGridOrbitals* grad_;
    LocGridOrbitals* pcgrad_;
    LocGridOrbitals* sdir_;

    static Timer line_min_tm_;
    static Timer nl_update_tm_;
    static Timer comp_res_tm_;
    static Timer update_states_tm_;

    void update_states(LocGridOrbitals& orbitals, LocGridOrbitals& res,
        LocGridOrbitals& work_orbitals, const double precond_factor);

public:
    GrassmanLineMinimization(Hamiltonian* hamiltonian,
        ProjectedMatricesInterface* proj_matrices, MGmol* mgmol_strategy,
        Ions& ions, std::ostream& os)
        : hamiltonian_(hamiltonian),
          proj_matrices_(proj_matrices),
          mgmol_strategy_(mgmol_strategy),
          os_(os)
    {
        ptr2ions_ = &ions;
        ;
        conjugate_ = false;
    }

    ~GrassmanLineMinimization()
    {
        delete new_grad_;
        delete new_pcgrad_;
        delete grad_;
        delete pcgrad_;
        delete sdir_;
    }

    void setup(LocGridOrbitals&){};

    int update(LocGridOrbitals& orbitals, Ions& ions,
        const double precond_factor, const bool orthof,
        LocGridOrbitals& work_orbitals, const bool accelerate,
        const bool print_res, const double atol);

    virtual void conjugate()                                  = 0;
    virtual double computeStepSize(LocGridOrbitals& orbitals) = 0;
    virtual void parallelTransportUpdate(
        const double lambda, LocGridOrbitals& orbitals)
        = 0;
    static void printTimers(ostream& os);
};
#endif
