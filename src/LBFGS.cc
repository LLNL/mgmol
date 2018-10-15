// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "LBFGS.h"
#include "Control.h"
#include "ProjectedMatrices.h"
#include "Mesh.h"
#include "MasksSet.h"
#include "Electrostatic.h"
#include "KBPsiMatrix.h"
#include "MGmol.h"


LBFGS::LBFGS(LocGridOrbitals** orbitals,
             Ions& ions,
             Rho& rho,
             ConstraintSet& constraints,
             LocalizationRegions& lrs,
             ClusterOrbitals* local_cluster,
             MasksSet& masks,
             MasksSet& corrmasks,
             Electrostatic& electrostat,
             const double dt,
             MGmol& strategy):
                IonicAlgorithm(orbitals,
                               ions,
                               rho,
                               constraints,
                               lrs,
                               masks,
                               strategy),
    orbitals_(orbitals),
    ions_(ions),
    rho_(rho),
    lrs_(lrs),
    local_cluster_(local_cluster),
    masks_(masks),
    corrmasks_(corrmasks),
    electrostat_(electrostat),
    ref_lrs_(lrs),
    mgmol_strategy_(strategy)
{
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();
    Control& ct = *(Control::instance());

    etot_i_[0]=10000.;
    etot_i_[1]=10000.;
    etot_i_[2]=10000.;

    stepper_
        = new LBFGS_IonicStepper(dt, atmove_, tau0_, taup_, fion_, gid_,
                                 20, &etot_i_[0]);
    registerStepper(stepper_);

    ref_masks_    =new MasksSet(lrs,false,ct.getMGlevels());
    ref_corrmasks_=new MasksSet(lrs,true,0);
    
    vh_init_  =new pb::GridFunc<POTDTYPE> (electrostat_.getVh());
    
    ref_orbitals_
        =new LocGridOrbitals(mygrid, 
                             mymesh->subdivx(),
                             ct.numst,
                             ct.bc,
                             (*orbitals_)->getProjMatrices(),
                             &ref_lrs_,
                             ref_masks_,ref_corrmasks_,local_cluster_);
    
    ref_orbitals_->setup(&ref_lrs_);
    ref_orbitals_->assign(**orbitals_);
}

LBFGS::~LBFGS()
{
    delete vh_init_;
    delete ref_masks_;
    delete ref_corrmasks_;
    delete stepper_;
    delete ref_orbitals_;
}

int LBFGS::quenchElectrons(const int itmax, double& etot)
{
    etot_i_[0] = etot_i_[1];
    etot_i_[1] = etot_i_[2];
    int ret=mgmol_strategy_.quench(*orbitals_, ions_, itmax, 0, etot_i_[2]);
    
    etot=etot_i_[2];
    return ret;
}

void LBFGS::updateRefs()
{
    Control& ct = *(Control::instance());
    if( !stepper_->check_last_step_accepted() )
    {
        if( onpe0 )
            (*MPIdata::sout)<<"lbfgs: Reset orbitals to reference orbitals "<<endl;
        (*orbitals_)->assign(*ref_orbitals_);
        electrostat_.setupInitialVh(*vh_init_);
    }
    else
    {
        // update references
        if( ct.lr_update )
        {
            ref_lrs_    = lrs_;
            *ref_masks_ = masks_;
            *ref_corrmasks_ = corrmasks_;
        }
        ref_orbitals_->reset(ref_masks_, ref_corrmasks_, &lrs_);
        if( onpe0 )
             (*MPIdata::sout)<<"lbfgs: Update reference orbitals "<<endl;
        ref_orbitals_->assign(**orbitals_);
        vh_init_->assign(electrostat_.getVh(),'d');
    }
}

void LBFGS::setQuenchTol()const
{
    Control& ct = *(Control::instance());
    ct.conv_tol=0.1*stepper_->etol();
    if( onpe0 )
        (*MPIdata::sout)<<setprecision(12)<<fixed
            <<"lbfgs: Set SC convergence criterion to "<<ct.conv_tol<<endl;
}

void LBFGS::updatePotAndMasks()
{
    if( !stepper_->check_last_step_accepted() )
        electrostat_.setupInitialVh(*vh_init_);

    IonicAlgorithm::updatePotAndMasks();
}

bool LBFGS::lbfgsLastStepNotAccepted()const
{
    return !stepper_->check_last_step_accepted();
} 
