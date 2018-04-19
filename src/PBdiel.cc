// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include <iostream>
#include <iomanip>
using namespace std;

#include "PBdiel.h"
#include "MGmol_blas1.h"


template class PBdiel<pb::PBh2<POTDTYPE> >;
//template class PBdiel<pb::PBh2<float> >;
template class PBdiel<pb::PBh4<POTDTYPE> >;
//template class PBdiel<pb::PBh4<float> >;
template class PBdiel<pb::PBh6<POTDTYPE> >;
//template class PBdiel<pb::PBh6<float> >;
template class PBdiel<pb::PBh8<POTDTYPE> >;
//template class PBdiel<pb::PBh8<float> >;
template class PBdiel<pb::PBh4M<POTDTYPE> >;
//template class PBdiel<pb::PBh4M<float> >;
template class PBdiel<pb::PBh4MP<POTDTYPE> >;
//template class PBdiel<pb::PBh4MP<float> >;

template <class T>
void PBdiel<T>::solve(const pb::GridFunc<RHODTYPE>& rho, 
                      const pb::GridFunc<RHODTYPE>& rhoc)
{
    PoissonInterface::poisson_tm_.start();

    pb::GridFunc<RHODTYPE> work_rho(rho);
  
    // Keep in memory vh*rho before updating vh 
    const double vel=Poisson::grid_.vel();
    Poisson::Int_vhrho_old_ = vel * Poisson::vh_->gdot(rho);

    // Subtract compensating charges from rho 
    work_rho-=rhoc;
    
    assert( rhod_!=NULL );

    
    /* Check for uniform precision before calling poisson_solver.
     * Downgrade or upgrade rhs (work_rho) and rhod_ to have precision of solution (vh_).
     * Note that this could be done at the beginning of this function, but several 
     * operations involving rho might be done in lower precision (depending on POTDTYPE), 
     * which could affect accuracy. For now, we delay the switch until just before the 
     * solve call. 
    */
//    if(sizeof(POTDTYPE) != sizeof(RHODTYPE))
//    {
       /* solve with POTDTYPE precision */
       pb::GridFunc<POTDTYPE> rhs(work_rho);
       pb::GridFunc<POTDTYPE> rhod(*rhod_);
       poisson_solver_->solve(*Poisson::vh_, rhs, rhod, *Poisson::vepsilon_);
//    }
//    else
//    {
//       poisson_solver_->solve(*Poisson::vh_, work_rho, *rhod_, *Poisson::vepsilon_);
//    }
   
    double residual_reduction = poisson_solver_->getResidualReduction();
    double final_residual     = poisson_solver_->getFinalResidual();

    if( onpe0 )
       (*MPIdata::sout)<<setprecision(2)<<scientific
           <<"PBdiel: residual reduction = "<<residual_reduction
           <<", final residual = "<<final_residual<<endl;

    Poisson::Int_vhrho_     = vel * Poisson::vh_->gdot(rho);
    Poisson::Int_vhrhoc_    = vel * Poisson::vh_->gdot(rhoc);

    PoissonInterface::poisson_tm_.stop();
}

template <class T>
void PBdiel<T>::set_rhod(pb::GridFunc<RHODTYPE>* rhod)
{
    //(*MPIdata::sout)<<"set_rhod"<<endl;
    assert( rhod!=NULL );
    rhod_=rhod;
}
