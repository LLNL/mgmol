// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "MGmol.h"
#include "Ions.h"
#include "Hamiltonian.h"
#include "AndersonMix.h"
#include "Potentials.h"
#include "Control.h"
#include "ABPG.h"
#include "ProjectedMatricesInterface.h"

bool ABPG::pbset_ = false;

Timer ABPG::abpg_tm_("abpg_line_min");
Timer ABPG::abpg_nl_update_tm_("abpg_nl_update");
Timer ABPG::comp_res_tm_("abpg_comp_residuals_st");
Timer ABPG::update_states_tm_("abpg_update_states");

void ABPG::setup(LocGridOrbitals& orbitals)
{
    Control& ct = *(Control::instance());
    
    if( ct.wf_dyn==1 )// use Anderson extrapolation
    {
        wf_mix_=new AndersonMix<LocGridOrbitals>(ct.wf_m, ct.betaAnderson, 
                                      orbitals, (ct.orbital_type==2) );    
    }
}

//
// Performs a single wave functions update step.
//
// orthof=true: wants orthonormalized updated wave functions
int ABPG::update(LocGridOrbitals& orbitals,
        Ions& ions,
        const double precond_factor,
        const bool orthof,
        LocGridOrbitals& work_orbitals,
        const bool accelerate,
        const bool print_res,
        const double atol)
{
    abpg_nl_update_tm_.start();
    
    Control& ct = *(Control::instance());
    if( onpe0 && ct.verbose>2 ) os_<<"ABPG::update()..."<<endl;
    
    // temporary LocGridOrbitals to hold residual
    LocGridOrbitals res(orbitals,false);

    const bool check_res=(atol>0.);
    double normRes=mgmol_strategy_->computeResidual(orbitals,work_orbitals,res,
       (print_res ||check_res),
       check_res );
    if( normRes<atol && check_res )
    {
        abpg_nl_update_tm_.stop();
        return 0;
    }

    // update wavefunctions
    update_states(orbitals, res, work_orbitals, precond_factor, accelerate);

    bool flag_ortho=false;
    if( ct.orbital_type==0 || orthof )
    {
        if( ct.isLocMode() )
        {
            //orbitals.ortho_norm_local();
            orbitals.normalize();
        }
        else
        {
            orbitals.orthonormalize();
            flag_ortho=true;
        }
        if( ct.wf_dyn ) // restart mixing after (ortho-)normalization
        {
            if(wf_mix_!=0)wf_mix_->restart();
        }
    }
    else
    if( ct.orbital_type!=2 )
    {
        orbitals.normalize();
    }
    
    // if orthonorm() not called, recompute overlap
    if( !flag_ortho && ct.orbital_type!=2 )
    {
        orbitals.computeGramAndInvS();
    }
    
#if 0
    int msize=ct.numst;
    if( ct.loc_mode ){
        // work matrices
        dist_matrix::DistMatrix<DISTMATDTYPE>  u_dis("U", msize, msize);
        u_dis=proj_matrices_->matS();
        dist_matrix::DistMatrix<DISTMATDTYPE>  w_dis("W", msize, msize);
        w_dis.identity();
        u_dis.axpy(-1.,w_dis);
        const double inorm=u_dis.norm('i');
        if( onpe0 && inorm>1.e-8 )
            os_<<setprecision(2)<<scientific
                <<"Largest non diagonal element of S = "<<inorm<<endl;
    }
#endif

    abpg_nl_update_tm_.stop();
    
    return 1;
} 

//////////////////////////////////////////////////////////////////////////////
// Update orbitals using MG preconditioning and Anderson extrapolation
void ABPG::update_states(LocGridOrbitals& orbitals, 
                         LocGridOrbitals& res,
                         LocGridOrbitals& work_orbitals,
                         const double precond_factor,
                         const bool accelerate)
{
    assert( orbitals.getIterativeIndex()>=0 );
    
    update_states_tm_.start();

    //if(onpe0) os_<<"Update wave functions"<<endl;

    Control& ct = *(Control::instance());

    if( (ct.getPrecondType()%10)==0 && ct.getMGlevels()>=0 )
    {
        // PRECONDITIONING 
        // compute the preconditioned steepest descent direction
        // -> res
        mgmol_strategy_->precond_mg(res);
    }
    
    //non-energy compatible spread penalties are added after 
    //preconditioning is applied
    if( ct.isSpreadFunctionalActive() && !ct.isSpreadFunctionalEnergy() )
        mgmol_strategy_->addResidualSpreadPenalty(orbitals,res);

    // apply AOMM projector to have a gradient orthogonal to kernel functions
    if( ct.use_kernel_functions )
    {
        mgmol_strategy_->projectOutKernel(res);
        
        // just reapply mask for now...
        res.applyMask();
    }

    if( ct.project_out_psd )
    {
        if(onpe0)
            os_<<"Project out preconditioned gradient"<<endl;
        orbitals.projectOut(res);
    }

    const double alpha=0.5*precond_factor;

    // Update wavefuntion
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        os_<<"Update states"<<endl;    
#else
    if( onpe0 && ct.verbose>2 )
        os_<<"Update states"<<endl;    
#endif
    if( accelerate )
    {
        res.scal(alpha);
        // Extrapolation scheme
        assert( wf_mix_!=0 );
        wf_mix_->update(res,work_orbitals);
    }
    else
    {
        // Preconditioned Power Method
        orbitals.axpy(alpha,res);
        
        if( ct.orbital_type==2 )
            orbitals.orthonormalizeLoewdin(false);
    }
    orbitals.incrementIterativeIndex();

    update_states_tm_.stop();

    assert( orbitals.getIterativeIndex()>=0 );
}

void ABPG::printTimers(ostream& os)
{
   abpg_tm_.print(os);
   abpg_nl_update_tm_.print(os);
   comp_res_tm_.print(os);
   update_states_tm_.print(os);      
}
