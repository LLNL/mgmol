// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "OrbitalsExtrapolationOrder3.h"
#include "LocGridOrbitals.h"
#include "MatricesBlacsContext.h"
#include "DistMatrixTools.h"
#include "ProjectedMatricesInterface.h"

#define EXTRAPOLATE_H 1

void OrbitalsExtrapolationOrder3::extrapolate_orbitals(LocGridOrbitals** orbitals, LocGridOrbitals* new_orbitals)
{
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext* bc = mbc. bcxt();

    Control& ct = *(Control::instance());

    new_orbitals->assign(**orbitals);

#if EXTRAPOLATE_H
    ProjectedMatricesInterface* proj_matrices=(*orbitals)->getProjMatrices();
#endif
    
    // do the extrapolation if previous orbitals exist (not at first step)
    
    if( orbitals_minus1_!=0 )
    {
        LocGridOrbitals tmp_orbitals_minus1(*new_orbitals, false);
        
        if( ct.verbose>1 && onpe0 )
            (*MPIdata::sout)<<"Extrapolate orbitals using 3rd order scheme..."<<endl;
        // align orbitals_minus1 with new_orbitals
        if( ct.it_algo_type>1 )
        {
            dist_matrix::DistMatrix<DISTMATDTYPE> matQ("Q",*bc, ct.numst, ct.numst);
            dist_matrix::DistMatrix<DISTMATDTYPE> yyt("yyt",*bc, ct.numst, ct.numst);

            // alignement
            orbitals_minus1_->computeGram(*new_orbitals, matQ);                
            getProcrustesTransform(matQ, yyt);
            orbitals_minus1_->multiply_by_matrix(matQ);

            // compute delta Phi
            tmp_orbitals_minus1.assign(*orbitals_minus1_);
            tmp_orbitals_minus1.axpy(-1.,*new_orbitals);
            tmp_orbitals_minus1.multiply_by_matrix(yyt);

#if EXTRAPOLATE_H
            proj_matrices->updateHminus1tmp(matQ,yyt,(*MPIdata::sout));
#endif
            
            if( orbitals_minus2_!=0 ){
#if 0
                LocGridOrbitals tmp(*orbitals_minus2_);
                tmp.axpy(-1.,*new_orbitals);
                tmp.computeGram(matQ);
                double normQ=matQ.trace();
                if( onpe0 )
                    (*MPIdata::sout)<<"||Phi_old-Phi_new|| before Procrustes = "<<normQ<<endl;
#endif                
                // alignement
                orbitals_minus2_->computeGram(*new_orbitals, matQ);                
                getProcrustesTransform(matQ, yyt);
                orbitals_minus2_->multiply_by_matrix(matQ);

                // compute delta Phi
                orbitals_minus2_->axpy(-1.,*orbitals_minus1_);
                orbitals_minus2_->multiply_by_matrix(yyt);
#if EXTRAPOLATE_H
                proj_matrices->updateHminus2(matQ,yyt, (*MPIdata::sout));
#endif

#if 0        
                tmp.assign(*orbitals_minus2_);
                tmp.computeGram(matQ);
                normQ=matQ.trace();
                if( onpe0 )
                    (*MPIdata::sout)<<"||Phi_old-Phi_new|| after Procrustes = "<<normQ<<endl;
#endif
            }
                      
#if EXTRAPOLATE_H
            if( ct.verbose>2 && onpe0 )
                (*MPIdata::sout)<<"Extrapolate H..."<<endl;
            proj_matrices->extrapolateHorder3();
#endif
        }else{
            tmp_orbitals_minus1.assign(*orbitals_minus1_);
            if( orbitals_minus2_!=0 )
                orbitals_minus2_->axpy(-1.,*orbitals_minus1_);
            if( ct.verbose>1 && onpe0 )
                (*MPIdata::sout)<<"Compute tmp_orbitals_minus1..."<<endl;
            tmp_orbitals_minus1.axpy(-1.,*new_orbitals);
        }
        
        if( orbitals_minus2_!=0 ){
            new_orbitals->axpy(-2.,tmp_orbitals_minus1);
            new_orbitals->axpy(1.,*orbitals_minus2_);
            
            delete orbitals_minus2_;
        }else{
            if( ct.verbose>1 && onpe0 )
                (*MPIdata::sout)<<"Extrapolate orbitals using 2nd order scheme only for this step..."<<endl;
            new_orbitals->axpy(-1.,tmp_orbitals_minus1);
        }

        orbitals_minus2_=orbitals_minus1_;
    }

    orbitals_minus1_=*orbitals;

    if( ct.it_algo_type>1 ){
#if EXTRAPOLATE_H
        proj_matrices->updateHminus2();
#endif
    }
    
    
    *orbitals=new_orbitals;
#if EXTRAPOLATE_H
    if( ct.it_algo_type>1 ){
        proj_matrices->saveH();
    }
#endif

    (*orbitals)->incrementIterativeIndex();

    if( ct.isLocMode() ){
        (*orbitals)->normalize();
        (*orbitals)->applyMask();
    }else{
        (*orbitals)->orthonormalizeLoewdin();
    }
    
}
