// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef _ABPG_H_
#define _ABPG_H_

#include "Mixing.h"
#include "OrbitalsStepper.h"
#include "LocGridOrbitals.h"

#include <iostream>

class ProjectedMatricesInterface;
class MGmol;
class Hamiltonian;
class Energy;

class ABPG:public OrbitalsStepper
{
    Hamiltonian* hamiltonian_;
    ProjectedMatricesInterface* proj_matrices_;
    Energy* energy_;
    MGmol* mgmol_strategy_;
    
    std::ostream& os_;
    
    Mixing<LocGridOrbitals>* wf_mix_;
    static bool pbset_;

    static Timer abpg_tm_;
    static Timer abpg_nl_update_tm_;
    static Timer comp_res_tm_;
    static Timer update_states_tm_;

    void update_states(LocGridOrbitals& orbitals, 
                       LocGridOrbitals& res,
                       LocGridOrbitals& work_orbitals,
                       const double precond_factor,
                       const bool accelerate);

public:

    ABPG(Hamiltonian* hamiltonian,
         ProjectedMatricesInterface* proj_matrices,
         MGmol* mgmol_strategy,
         std::ostream& os)
        : hamiltonian_(hamiltonian),
          proj_matrices_(proj_matrices),
          mgmol_strategy_(mgmol_strategy),
          os_(os),
          wf_mix_(0)
    {
    }
    
    ~ABPG()
    {
        if( wf_mix_!=0 )
        {
            delete wf_mix_;
            wf_mix_=0;
        }
    }
    
    void setup(LocGridOrbitals&);
    
    int update(LocGridOrbitals& orbitals,
        Ions& ions,
        const double precond_factor,
        const bool orthof,
        LocGridOrbitals& work_orbitals,
        const bool accelerate,
        const bool print_res,
        const double atol);

    static void printTimers(std::ostream& os);
};

#endif  
