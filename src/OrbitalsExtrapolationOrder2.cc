// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "OrbitalsExtrapolationOrder2.h"
#include "Control.h"
#include "LocGridOrbitals.h"
#include "DistMatrixTools.h"
#include "ProjectedMatrices.h"

#define EXTRAPOLATE_H 1

void OrbitalsExtrapolationOrder2::extrapolate_orbitals(LocGridOrbitals** orbitals, 
                                                       LocGridOrbitals* new_orbitals)
{
    Control& ct = *(Control::instance());
    
#if EXTRAPOLATE_H
    ProjectedMatricesInterface* proj_matrices=(*orbitals)->getProjMatrices();
    ProjectedMatrices* projmat=0;
    if( ct.it_algo_type>1 )
    {
        projmat = dynamic_cast<ProjectedMatrices*>(proj_matrices);
        assert( projmat );
    }
#endif

    new_orbitals->assign(**orbitals);

    extrapolated_H_=false;
    // do the extrapolation if previous orbitals exist (not at first step)
    {
        if( orbitals_minus1_!=0 )
        {
            if( ct.verbose>1 && onpe0 )
                (*MPIdata::sout)<<"Extrapolate orbitals order 2..."<<endl;
            
            // align orbitals_minus1_ with new_orbitals
            if( ct.it_algo_type>1 )
            {
                dist_matrix::DistMatrix<DISTMATDTYPE> matQ("Q", ct.numst, ct.numst);
#if 0
                LocGridOrbitals tmp(*orbitals_minus1_);
                tmp.axpy(-1.,*new_orbitals);
                tmp.computeGram(matQ);
                double normQ=matQ.trace();
                if( onpe0 )
                    (*MPIdata::sout)<<"||Phi_old-Phi_new|| before Procrustes = "<<normQ<<endl;
#endif                
                orbitals_minus1_->computeGram(*new_orbitals, matQ);
                
                dist_matrix::DistMatrix<DISTMATDTYPE> yyt("yyt", ct.numst, ct.numst);
                getProcrustesTransform(matQ, yyt);
                
                orbitals_minus1_->multiply_by_matrix(matQ);
                orbitals_minus1_->axpy(-1.,*new_orbitals);
                orbitals_minus1_->multiply_by_matrix(yyt);
                
#if EXTRAPOLATE_H
                if( ct.verbose>2 && onpe0 )
                    (*MPIdata::sout)<<"Extrapolate H..."<<endl;
                projmat->extrapolateHorder2(matQ,yyt,(*MPIdata::sout));
                extrapolated_H_=true;
#endif
            
#if 0          
                tmp.assign(*orbitals_minus1_);
                //tmp.axpy(-1.,*new_orbitals);
                tmp.computeGram(matQ);
                normQ=matQ.trace();
                if( onpe0 )
                    (*MPIdata::sout)<<"||Phi_old-Phi_new|| after Procrustes = "<<normQ<<endl;
#endif                
            }
            else
            { // ct.it_algo_type==0/1
            
                new_orbitals->scal(2.);
            }
            new_orbitals->axpy(-1.,*orbitals_minus1_);
            
            delete orbitals_minus1_;
        }
    
        // save data for next extrapolation
        if( ct.verbose>1 && onpe0 )
            (*MPIdata::sout)<<"Set orbitals_minus1_ to values of orbitals"<<endl;
        orbitals_minus1_=*orbitals;

        if( ct.it_algo_type>1 )
        {
#if EXTRAPOLATE_H
            ProjectedMatrices* projmat =
                dynamic_cast<ProjectedMatrices*>(proj_matrices);
            assert( projmat );

            projmat->updateHminus1();
#endif
        }
    }
    
    *orbitals=new_orbitals;
#if EXTRAPOLATE_H
    if( ct.it_algo_type>1 )
    {
        projmat->saveH();
    }
#endif
    
    (*orbitals)->incrementIterativeIndex();

    if( ct.isLocMode() )
    {
        (*orbitals)->normalize();
        (*orbitals)->applyMask();
    }
    else
    {
        // DM (if not recomputed from scratch)
        // is consistant with orthonormal set of orbitals... 
        (*orbitals)->orthonormalizeLoewdin();
    }
}
