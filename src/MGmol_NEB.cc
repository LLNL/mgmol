// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "MGmol.h"
#include "LBFGS.h"
#include "FIRE.h"
#include "Energy.h"
#include "DFTsolver.h"
#include "ProjectedMatricesInterface.h"

void MGmol::sebprintForces()
{
   ions_->printForces(os_);
}
void MGmol::sebprintPositions()
{
    ions_->printPositions(os_);
}

void MGmol::geomOptimSetup()
{
    Control& ct = *(Control::instance());

    switch(ct.atoms_dyn)
    {
        case 6:
            geom_optimizer_=new LBFGS(&current_orbitals_, *ions_, *rho_, *constraints_, 
                *lrs_, local_cluster_, *currentMasks_, *corrMasks_, *electrostat_, ct.dt,
                *this);
            break;
        
        case 7:
            geom_optimizer_=new FIRE(&current_orbitals_, *ions_, *rho_, *constraints_,
                *lrs_, *currentMasks_, *electrostat_, ct.dt,
                *this);
            break;
            
        default:
            (*MPIdata::serr)<<"geomOptimSetup(): option "<<ct.atoms_dyn
                            <<" is an invalid method"<<endl;
            return;
    }  
    DFTsolver::resetItCount();

    geom_optimizer_->init(h5f_file_);

    // additional quench to compensate random start
    if(ct.restart_info<3)
    {
        double eks=0.;
        geom_optimizer_->quenchElectrons(ct.max_electronic_steps, eks);
    }
    else
    {
        DFTsolver::setItCountLarge();
    }    
}

void MGmol::geomOptimQuench()
{
    Control& ct = *(Control::instance());

    double eks=0.;
    geom_optimizer_->quenchElectrons(ct.max_electronic_steps, eks);

    // Get the total energy
    double ts=0.5*proj_matrices_->computeEntropy(); // in [Ha]
    total_energy_ = energy_->evaluateTotal(ts, proj_matrices_, *current_orbitals_, 2, os_);
}

void MGmol::geomOptimComputeForces()
{
    geom_optimizer_->computeForces();
}

void MGmol::geomOptimSetForces(const vector<vector<double> >& f)
{
    geom_optimizer_->setForces(f);
}

void MGmol::geomOptimDumpRestart()
{
    geom_optimizer_->dumpRestart();
}

int MGmol::geomOptimRun1Step()
{
    int conv=geom_optimizer_->run1step();
    geom_optimizer_->updatePotAndMasks();
    // Write down positions and displacements
    ions_->printPositions(os_);    
    return conv;
}

short MGmol::geomOptimCheckTolForces(const double tol_force)
{
    return geom_optimizer_->checkTolForces(tol_force);
}
