// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Control.h"
#include "NonOrthoDMStrategy.h"
#include "ProjectedMatricesInterface.h"
#include "LocGridOrbitals.h"

NonOrthoDMStrategy::NonOrthoDMStrategy(
    LocGridOrbitals* orbitals,
    ProjectedMatricesInterface* proj_matrices,
    const double mix):
        orbitals_(orbitals),
        proj_matrices_(proj_matrices),
        mix_(mix)
{
}

void NonOrthoDMStrategy::initialize()
{
    Control& ct = *(Control::instance());
    
    if( onpe0 && ct.verbose>2 )
    {
        (*MPIdata::sout)<<"NonOrthoDMStrategy::initialize()..."<<endl;
    }
    proj_matrices_->updateDM(orbitals_->getIterativeIndex());        
}

int NonOrthoDMStrategy::update()
{
    assert( proj_matrices_!=0 );

    Control& ct = *(Control::instance());
    if( onpe0 && ct.verbose>2 )
    {
        (*MPIdata::sout)<<"NonOrthoDMStrategy::update() with mixing = "
                        <<mix_<<endl;
    }

    // save old density matrix
    if( mix_<1. )proj_matrices_->saveDM();
      
    // compute new density matrix
    proj_matrices_->updateDM(orbitals_->getIterativeIndex()); 
    
    if( mix_<1. )
    {
        proj_matrices_->updateDMwithRelax(mix_,orbitals_->getIterativeIndex());
    }

    if( ct.verbose>2 )
    {
        double dd=proj_matrices_->getNel();
        if( onpe0 )
            (*MPIdata::sout)<<setprecision(8)
                            <<"test NonOrthoDMStrategy::update(): Nel = "<<dd<<endl;
    }
    
    return 0; // success
}

void NonOrthoDMStrategy::stripDM()
{
    if( mix_<1. )
        proj_matrices_->stripDM();
}

void NonOrthoDMStrategy::dressDM()
{
    if( mix_<1. )
        proj_matrices_->dressupDM();
}

void NonOrthoDMStrategy::reset()
{
    proj_matrices_->resetDM();
}
