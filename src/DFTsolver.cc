// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DFTsolver.h"
#include "Control.h"
#include "MGmol.h"
#include "ABPG.h"
#include "Hamiltonian.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "ProjectedMatricesInterface.h"
#include "GrassmanCGSparse.h"
#include "GrassmanCG.h"
#include "Rho.h"
#include "DMStrategy.h"
#include "Ions.h"

Timer DFTsolver::solve_tm_("solve");
int DFTsolver::it_scf_=0;

DFTsolver::DFTsolver(Hamiltonian* hamiltonian,
                     ProjectedMatricesInterface* proj_matrices,
                     Energy* energy,
                     Electrostatic* electrostat,
                     MGmol* mgmol_strategy,
                     Ions& ions,
                     Rho* rho,
                     DMStrategy* dm_strategy,
                     std::ostream& os):
    hamiltonian_(hamiltonian),
    proj_matrices_(proj_matrices),
    energy_(energy),
    electrostat_(electrostat),
    mgmol_strategy_(mgmol_strategy),
    ions_(ions),
    rho_(rho),
    dm_strategy_(dm_strategy),
    os_(os)
{
    Control& ct(*(Control::instance()));
    
    switch ( ct.it_algo_type )
    {
        case 0:
        {
            orbitals_stepper_=new ABPG(hamiltonian_,
                          proj_matrices_, 
                          mgmol_strategy,
                          os_);
            break;
        }

        case 1:
        {
            if(ct.short_sighted)
            {
                orbitals_stepper_=new GrassmanCGSparse(hamiltonian_,
                        proj_matrices_, 
                        mgmol_strategy,
                        ions,
                        os_);
            }
            else
            {
                orbitals_stepper_=new GrassmanCG(hamiltonian_,
                        proj_matrices_, 
                        mgmol_strategy,
                        ions,
                        os_);
            }            
            break;
        }
        
        default:
            cerr<<"DFTsolver: Undefined iterative electronic structure solver!!!"<<endl;
    }
    
    accelerate_=false;
    //if(ct.restart_info>2 && ct.atoms_dyn == 2 && ct.wf_dyn==1 )accelerate_=true;
}

DFTsolver::~DFTsolver()
{
    delete orbitals_stepper_;
}

void DFTsolver::printEnergy(const short step)const
{   
    if ( onpe0 )
    {
        os_<<setprecision(12)<<fixed
           <<"%%    "<<step<<" SC ENERGY = "<<eks_history_[0];
        if( step>0 )
        {
            if( testUpdatePot() )
            {
                os_<<setprecision(2)<<scientific<<", delta Eks = "<<de_<<endl;
            }
            else
            {
                os_<<setprecision(2)<<scientific<<", delta Eig. sum = "<<deig_<<endl;
                os_<<setprecision(12)<<fixed
                    <<"Sum eigenvalues = "<<sum_eig_[0]<<endl;
            }
        }
        else
        {
            os_<<endl;
        }
    }
}

bool DFTsolver::checkPrintResidual(const short step)const
{
    Control& ct(*(Control::instance()));
    return (ct.iprint_residual>0) ? !(step%ct.iprint_residual) : false;
}

void DFTsolver::dielON()
{
    Potentials& pot(hamiltonian_->potential());    
    bool isON=pot.diel();
    if( !isON )return; // continuum solvent is OFF
    
    static bool pbset=false;
    if( pbset )return; // continuum solvent already set

    Control& ct(*(Control::instance()));
    const int diel_delay = (ct.restart_info<3) ? 10 : 1;
    if( it_scf_<diel_delay )
    {
        isON=false;
    }
    
    // turn ON continuum solvent
    if( isON && deig2_ < 2.e-2*ct.numst )
    {
        electrostat_->setupPB(ct.rho0,ct.drho0,pot);
        electrostat_->setup(ct.vh_its);
        pbset=true;
    }

    if( pot.diel() && !pbset )
    if( onpe0 && ct.verbose>1 )
        os_<<" Solvation turned off for this step"<<endl;
}

bool DFTsolver::testUpdatePot()const
{
    Control& ct(*(Control::instance()));
    return ( it_scf_>ct.max_changes_pot );
}

bool DFTsolver::checkConvPot()const
{
    Control& ct(*(Control::instance()));
    Potentials& pot(hamiltonian_->potential());
    if ( (fabs(pot.scf_dvrho()/ions_.getNumIons()) < ct.conv_tol) 
      && (it_scf_>(2+ct.max_changes_pot) ) ) 
        return true;

    return false;
}

// returns:
// 0 if converged, 
// 1 if not converged yet,
// -1 if not reaching minimum convergence
// -2 if failing to converge
int DFTsolver::checkConvergenceEnergy(const short step, const short max_steps)
{
    Control& ct ( *(Control::instance()) );
    
    double dtol2=1000.;
    if( step>0 )
    {
        double deig_old=deig_;
        deig_    =fabs(sum_eig_[1]-sum_eig_[0]);
        deig2_   =max(deig_, deig_old);
        
        // energy variation during last iteration
        double de_old  =de_;
        de_     =fabs(eks_history_[0]-eks_history_[1]);
        de2_    =max(de_, de_old);

        if( testUpdatePot() ) dtol2=de2_;
        else                                   dtol2=deig2_;
    }
    
    // accelerate only if close enough to solution
    accelerate_ = ( accelerate_ ||
                 ( ct.wf_dyn && (deig2_ < 1.e-2*ct.numst) ) );
    
    // test if this iteration converged
    if( step>1 && dtol2<ct.conv_tol )return 0;

    if( step==max_steps && dtol2>ct.conv_tol_stop )return -1;

    // Test very bad behavior!!
    MGmol_MPI& mmpi ( *(MGmol_MPI::instance()) );
    Potentials& pot( hamiltonian_->potential() );
    if( pot.scf_dv() > 1. && step>50 && step>ct.max_changes_pot )
    {
        if( onpe0 )
        {
            os_<<endl
               <<"MGmol ERROR, electronic iterations are not converging at all: DV ="<<pot.scf_dv()
               <<endl<<endl;
            os_<<flush;
        }
        mmpi.barrier();
        return -2;
    }

    return 1;
}

double DFTsolver::evaluateEnergy(const LocGridOrbitals& orbitals, const bool print_flag)
{
    // save energy recent history
    eks_history_[1] = eks_history_[0];

    // Get the new total energy
    const double ts=0.5*proj_matrices_->computeEntropy(); // in [Ha]
    eks_history_[0] = energy_->evaluateTotal(ts, proj_matrices_, orbitals, print_flag, os_);
    
    sum_eig_[1] = sum_eig_[0];
    sum_eig_[0] = 2.*proj_matrices_->getEigSum(); // 2.*sum in [Ry]
    
    return eks_history_[0];
}

int DFTsolver::solve(LocGridOrbitals& orbitals,
         LocGridOrbitals& work_orbitals,
         Ions& ions, 
         const short max_steps,
         const short iprint,
         double& last_eks)
{
    solve_tm_.start();
    
    Control& ct(*(Control::instance()));

    eks_history_[0]=1.e9;
    eks_history_[1]=1.e9;
    
    sum_eig_[0]=1.e9;
    sum_eig_[1]=1.e9;
    
    int retval=1; // 0 -> converged, -1 -> problem, -2 -> ( de>conv_tol_stop )

    deig_  =1.e9;
    deig2_ =1.e9;
    de_    =1.e9;
    de2_   =1.e9;

#ifdef HAVE_ARPACK
    if( ct.precond_factor_computed )
    {
        const double ela=getLAeigen(0., 500, ions);
        if( fabs(ela)>1.e-16 )
            ct.precond_factor=1./ela;
        else
            return -1;
    } 
#else
    if( ct.precond_factor_computed )
    {
        os_<<"Needs ARPACK to compute Preconditioner factor"<<endl;
        return -1;
    }
#endif
    if( onpe0 )
    {
        os_<<"### DFTsolver ###"<<endl;
        if( ct.verbose>1 && ct.precond_factor_computed )
            os_<<"Preconditioning factor: "<<ct.precond_factor<<endl;
    }
   
    orbitals_stepper_->setup(orbitals);
   
    if( ct.resetVH() )
    {
        if( onpe0 )os_<<"DFTsolver: reset Hartree potential"<<endl;
        electrostat_->resetSolution();
    }
    
    // main electronic structure outer loop
    for(short step = 0;step <= max_steps;step++)
    {
        bool orthof=false;
        if( ct.orthof )
            orthof = ( ( ((step+1)%ct.orthof)==ct.max_changes_pot ) 
                       || ct.orthof==1 
                       || step==max_steps-1 );
        
        // turn on PB solver if necessary
        dielON();
        
        proj_matrices_->resetDotProductMatrices();
        
        // Generate new density
        rho_->update(orbitals);
        
        // Update potential
        if( testUpdatePot() ) mgmol_strategy_->update_pot(ions);
        
        mgmol_strategy_->updateHmatrix(orbitals,ions);

        // theta = invB * Hij
        // (to be used for energy and gradient computation)
        proj_matrices_->updateThetaAndHB();

        if( step == max_steps ) break;
        
        // Output the eigenvalues and occupations
        bool flag=false;
        if( iprint )flag=( !(step%iprint) && step<max_steps-1 );

        if( flag ) mgmol_strategy_->printEigAndOcc();

        last_eks = evaluateEnergy(orbitals,flag);

        // test for convergence
        retval=checkConvergenceEnergy(step, max_steps);
        
        // terminate if convergence problem
        if( retval<0 )return retval;

        printEnergy(step);

        // terminate early if convergence achieved
        if( retval==0 && !ct.checkResidual() )
        {
            if( onpe0 )os_<<endl<<endl<<" DFTsolver: convergence achieved for delta E..."<<endl;
            break;
        }
        
        bool print_res = checkPrintResidual(step);

        // strip dm from the overlap contribution
        // dm <- Ls**T * dm * Ls
        dm_strategy_->stripDM();

        // one step wave functions update
        // S and S^-1 should be up to date after that call
        const double restol = ct.checkResidual() ? ct.conv_tol : -1.;
        retval=orbitals_stepper_->update(orbitals, ions, ct.precond_factor, orthof, work_orbitals, 
                                         accelerate_, print_res, restol);

        // rebuild dm with new overlap matrix
        dm_strategy_->dressDM();

        if( retval==0 )
        {
            if( onpe0 )os_<<endl<<endl<<" DFTsolver: convergence achieved for residual..."<<endl;
            break;
        }


        // update ghost data now, so that it's ready to use later
        orbitals.setDataWithGhosts();
        orbitals.trade_boundaries();

        // compute B and its inverse for Mehrstellen, so that it's ready to use later
        // (depends on Phi only)
        orbitals.computeBAndInvB(*(hamiltonian_->lapOper()));
        
        if( ct.computeCondGramQuench() )
        {
            double condS=proj_matrices_->computeCond();
            if( onpe0 )
                os_<<setprecision(2)<<scientific<<"Condition Number of S: "
                   <<condS<<endl;
        }
        
        // updated Hij needed to compute new DM
        if( dm_strategy_->needH() )mgmol_strategy_->updateHmatrix(orbitals,ions);
        
        // compute new density matrix
        dm_strategy_->update();
        
    
        incInnerIt();

    } // end iterations
    
    if( iprint ) mgmol_strategy_->printEigAndOcc();

#ifdef HAVE_ARPACK
    if( ct.precond_factor_computed )
    {
        const double alpha = 1./getLAeigen(0., 500, ions);
        if( onpe0 ){
            os_<<"Quench electrons"<<endl;
            os_<<"A Posteriori Preconditioning factor: "
                <<alpha<<endl;
        }
    }
#endif
    solve_tm_.stop();

    return retval;
}

void DFTsolver::printTimers(ostream& os)
{
   solve_tm_.print(os);
}
