// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#include "HamiltonianMVPSolver.h"
#include "LocGridOrbitals.h"
#include "Ions.h"
#include "MatricesBlacsContext.h"
#include "Control.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesSparse.h"
#include "Rho.h"
#include "Energy.h"
#include "MGmol.h"
#include "KBPsiMatrixSparse.h"
#include "Electrostatic.h"
#include "Potentials.h"
#include "tools.h"
#include "DistMatrixWithSparseComponent.h"

#include <iomanip>
double evalEntropyMVP(ProjectedMatricesInterface* projmatrices,
                   const bool print_flag,
                   ostream& os);

template <class T1, class T2, class T3>
Timer HamiltonianMVPSolver<T1, T2, T3>::solve_tm_("HamiltonianMVPSolver::solve");
template <class T1, class T2, class T3>
Timer HamiltonianMVPSolver<T1, T2, T3>::target_tm_("HamiltonianMVPSolver::target");


//static int sparse_distmatrix_nb_partitions=128;

template <class T1, class T2, class T3>
HamiltonianMVPSolver<T1, T2, T3>::HamiltonianMVPSolver(MPI_Comm comm, ostream& os,
                           Ions& ions,
                           Rho* rho,
                           Energy* energy,
                           Electrostatic* electrostat,
                           MGmol* mgmol_strategy,
                           const int numst,
                           const double kbT,
                           const int nel,
                           const vector<vector<int> >& global_indexes,
                           const short n_inner_steps,
                           const T1& hinit,
                           const bool try_shorter_intervals):
    comm_(comm),
    os_(os),
    ions_(ions),
    n_inner_steps_(n_inner_steps),
    try_shorter_intervals_(try_shorter_intervals)
{
    assert( n_inner_steps>0 );
    
    rho_=rho;
    energy_=energy;
    electrostat_=electrostat;
    mgmol_strategy_=mgmol_strategy;
    
    numst_= numst;
    
    hmatrix_ = new T1(hinit);
    initial_hmatrix_ = new T1(hinit);
}

template <class T1, class T2, class T3>
HamiltonianMVPSolver<T1, T2, T3>::~HamiltonianMVPSolver()
{
    delete hmatrix_;
    delete initial_hmatrix_;
}

template <class T1, class T2, class T3>
void HamiltonianMVPSolver<T1, T2, T3>::reset()
{
    (*hmatrix_)=(*initial_hmatrix_);
}


// update density matrix in N x N space
template <class T1, class T2, class T3>
int HamiltonianMVPSolver<T1, T2, T3>::solve(LocGridOrbitals& orbitals)
{
    Control& ct = *(Control::instance());

    assert( numst_=(int)orbitals.numst() );
    assert( n_inner_steps_>0 );

    solve_tm_.start();
    
    if( onpe0 && ct.verbose>1 )
    {
        os_<<"----------------------------------------------------------------"<<endl;
        os_<<"Update DM functions using Hamiltonian MVP Solver..."<<endl;
        os_<<"----------------------------------------------------------------"<<endl;
    }
    
    //save initial matrix to enable reset
    (*initial_hmatrix_)=(*hmatrix_);

    KBPsiMatrixSparse kbpsi(0);
    kbpsi.setup(ions_, orbitals);
    
    T3 * projmatrices=dynamic_cast<T3 *>( orbitals.getProjMatrices() );

    int iterative_index=0;

    
    // save computed vh for a fair energy "comparison" with vh computed 
    // in close neigborhood
    const pb::GridFunc<POTDTYPE> vh_init(electrostat_->getVh());

    orbitals.setDataWithGhosts();
    
    // compute linear component of H
//    T2  h11("h11", numst_, numst_,comm_,
//                                 mgmol_strategy_->getRemoteTasksDistMatrix(),
//                                 sparse_distmatrix_nb_partitions);

    T2  h11("h11", numst_, comm_);

    kbpsi.computeAll(ions_, orbitals);
    
    kbpsi.computeHvnlMatrix(&kbpsi,ions_,h11);
    
    for(int inner_it=0;inner_it<n_inner_steps_;inner_it++)
    {
        if( onpe0 && ct.verbose>1 )
        {
            os_<<"---------------------------"<<endl;
            os_<<"Inner iteration "<<inner_it<<endl;
            os_<<"---------------------------"<<endl;
        }
    
        //
        // evaluate energy at origin
        //
        iterative_index++;
    
        projmatrices->assignH(*hmatrix_);
        projmatrices->setHB2H();

        // update DM and compute entropy
        projmatrices->updateDM(iterative_index); 
        double ts0 = evalEntropyMVP(projmatrices,true,os_);       
        // Update density
        rho_->update(orbitals);
        
        // Update potential
        mgmol_strategy_->update_pot(vh_init, ions_);

        energy_->saveVofRho();


        //compute new h11 for the current potential by adding local part to
        //nonlocal components (old local part reset to zero) 
        mgmol_strategy_->addHlocal2matrix(orbitals, orbitals, h11);

        projmatrices->assignH(h11);
        projmatrices->setHB2H();


        //compute energy at origin
        const int printE = ( ct.verbose>1 ) ? 1 : 0;
        double e0 = energy_->evaluateTotal(ts0, projmatrices, orbitals, printE, os_);
        
        //
        // compute energy at end for new H
        //
        
        T1 htarget(projmatrices->getH());

        iterative_index++;
        
        // update DM and compute entropy
        projmatrices->updateDM(iterative_index); 
        double ts1 = evalEntropyMVP(projmatrices,true,os_);
        // Update density
        rho_->update(orbitals);
        
        // Update potential
        mgmol_strategy_->update_pot(vh_init, ions_);

        energy_->saveVofRho();

        
        // update H and compute energy at midpoint
        mgmol_strategy_->addHlocal2matrix(orbitals, orbitals, h11);

        projmatrices->assignH(h11);
        projmatrices->setHB2H();


        //compute energy at end (beta=1.)
        double e1 = energy_->evaluateTotal(ts1, projmatrices, orbitals, printE, os_);


        //
        // evaluate energy at mid-point
        //
        T1 delta_h(htarget);
        delta_h -= *hmatrix_;
        
        h11 = *hmatrix_;
        h11.axpy(0.5,delta_h);

        projmatrices->assignH(h11);
        projmatrices->setHB2H();
    
        iterative_index++;
        
        //update DM and entropy
        projmatrices->updateDM(iterative_index); 
        double tsi = evalEntropyMVP(projmatrices,true,os_);
                
        // Update density
        rho_->update(orbitals);
        
        // Update potential
        mgmol_strategy_->update_pot(vh_init, ions_);

        energy_->saveVofRho();

        
        // update H with new potential
        mgmol_strategy_->addHlocal2matrix(orbitals, orbitals, h11);

        projmatrices->assignH(h11);
        projmatrices->setHB2H();


        //compute energy at midpoint
        double ei = energy_->evaluateTotal(tsi, projmatrices, orbitals, printE, os_);

        // line minimization
        double beta = minQuadPolynomialFrom3values(e0,e1,ei, (ct.verbose>2), os_);

        if( onpe0 && ct.verbose>1 )
        {
            os_<<setprecision(12);
            os_<<fixed<<"Inner iteration "<<inner_it<<", E0="<<e0<<", E(1/2)="<<ei<<", E1="<<e1;
            os_<<scientific<<" -> beta="<<beta;
            os_<<endl;
        }
        
        if( try_shorter_intervals_ )
        {
            double factor=0.5;
            while( beta<0. ) // try with a shorter interval if line search failed
            {
                if( onpe0 && ct.verbose>1 )
                {
                    os_<<"HMVP: Reduce interval by factor "<<factor<<" ..."<<endl;
                }
                ts1=tsi;
                e1 =ei;
            
                h11 = *hmatrix_;
                h11.axpy(0.5*factor,delta_h);
            
                projmatrices->assignH(h11);
                projmatrices->setHB2H();
            
                iterative_index++;
                
                // update DM and entropy
                projmatrices->updateDM(iterative_index); 
                tsi = evalEntropyMVP(projmatrices,true,os_);
                                
                // Update density
                rho_->update(orbitals);
                
                // Update potential
                mgmol_strategy_->update_pot(vh_init, ions_);
            
                energy_->saveVofRho();
               
                // update H            
                mgmol_strategy_->addHlocal2matrix(orbitals, orbitals, h11);
            
                projmatrices->assignH(h11);
                projmatrices->setHB2H();

            
                //compute energy at end (beta=1.)
                ei = energy_->evaluateTotal(tsi, projmatrices, orbitals, printE, os_);
            
                // line minimization
                beta = minQuadPolynomialFrom3values(e0,e1,ei, (ct.verbose>2), os_);
            
                if( onpe0 && ct.verbose>1 )
                {
                    os_<<setprecision(12);
                    os_<<fixed<<"Inner iteration "<<inner_it<<", E0="<<e0<<", E(1/2)="<<ei<<", E1="<<e1;
                    os_<<scientific<<" -> beta="<<beta;
                    os_<<endl;
                }
                
                beta*=factor;
                
                factor*=0.5;
            }
        }
        else
        {
            if( beta<0. )
            {
                if( onpe0 )os_<<"!!! HMVP iteration failed: beta<0 !!!"<<endl;
                projmatrices->assignH(*hmatrix_);
                projmatrices->setHB2H();
                
                return -1;
            }
        }

        hmatrix_->axpy(beta,delta_h);
        

    } // inner iterations
    
    projmatrices->assignH(*hmatrix_);
    projmatrices->setHB2H();

    iterative_index++;
    
    projmatrices->updateDM(iterative_index); 
    
    // Generate new density
    rho_->update(orbitals);
    
    if( onpe0 && ct.verbose>1 )
    {
        os_<<"----------------------------------------------------------------"<<endl;
        os_<<"End Hamiltonian MVP Solver..."<<endl;
        os_<<"----------------------------------------------------------------"<<endl;
    }
    solve_tm_.stop();
    
    return 0;
}

template <class T1, class T2, class T3>
void HamiltonianMVPSolver<T1, T2, T3>::printTimers(ostream& os)
{
    if( onpe0 )
    {
        os<<setprecision(2)<<fixed<<endl;
        solve_tm_.print(os);
        target_tm_.print(os);
    }
}

template class HamiltonianMVPSolver< dist_matrix::DistMatrix<DISTMATDTYPE>, dist_matrix::DistMatrixWithSparseComponent<DISTMATDTYPE>, ProjectedMatrices >;
template class HamiltonianMVPSolver< VariableSizeMatrix<sparserow>, VariableSizeMatrix<sparserow>, ProjectedMatricesSparse >;
//template int HamiltonianMVPSolver< dist_matrix::DistMatrix<DISTMATDTYPE>, dist_matrix::DistMatrixWithSparseComponent<DISTMATDTYPE> >::solve< ProjectedMatrices * >(LocGridOrbitals& orbitals);
