// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Control.h"
#include "MVP_DMStrategy.h"
#include "ProjectedMatricesInterface.h"
#include "LocGridOrbitals.h"
#include "MVPSolver.h"
#include "Ions.h"

#include <vector>
using namespace std;

MVP_DMStrategy::MVP_DMStrategy(
    MPI_Comm comm, 
    ostream& os, 
    Ions& ions,
    Rho* rho,
    Energy* energy,
    Electrostatic* electrostat,
    MGmol* mgmol_strategy,
    LocGridOrbitals* orbitals,
    ProjectedMatricesInterface* proj_matrices,
    const bool use_old_dm):
        comm_(comm),
        os_(os),
        ions_(ions),
        rho_(rho),
        energy_(energy),
        electrostat_(electrostat),
        global_indexes_(orbitals->getGlobalIndexes()),
        orbitals_(orbitals),
        proj_matrices_(proj_matrices),
        mgmol_strategy_(mgmol_strategy),
        use_old_dm_(use_old_dm)
{
    assert( electrostat_!=0 );
    assert( energy_!=0 );
}

void MVP_DMStrategy::initialize()
{
}

int MVP_DMStrategy::update()
{
    Control& ct = *(Control::instance());
    if( onpe0 && ct.verbose>2 )
    {
        (*MPIdata::sout)<<"MVP_DMStrategy::update()..."<<endl;
    }

    MVPSolver solver(comm_, os_,ions_,
             rho_, energy_, electrostat_,
             mgmol_strategy_,
             ct.numst,ct.occ_width,
             ct.getNel(),
             global_indexes_,
             ct.dm_inner_steps,
             use_old_dm_);
    
    return solver.solve(*orbitals_);
}

void MVP_DMStrategy::stripDM()
{
    if(use_old_dm_)proj_matrices_->stripDM();
}

void MVP_DMStrategy::dressDM()
{
    if(use_old_dm_)proj_matrices_->dressupDM();
}
