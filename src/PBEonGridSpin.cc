// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#include "PBEonGridSpin.h"

#include "Control.h"
#include "Mesh.h"
#include "Delh4.h"
#include "PBh4.h"

#include "Potentials.h"
#include "mputils.h"

PBEonGridSpin::PBEonGridSpin(Rho& rho, Potentials& pot):
    np_(rho.rho_[0].size()),
    rho_(rho),
    pot_(pot)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    myspin_=mmpi.myspin();
#ifdef USE_LIBXC
    int func_id = XC_GGA_X_PBE;
    if(xc_func_init(&xfunc_, func_id, XC_POLARIZED) != 0){
        cerr<<"Functional "<<func_id<<" not found"<<endl;
    }
    func_id = XC_GGA_C_PBE;
    if(xc_func_init(&cfunc_, func_id, XC_POLARIZED) != 0){
        cerr<<"Functional "<<func_id<<" not found"<<endl;
    }
    exc_.resize(np_*2);
    vsigma_.resize(np_*3);
#else
    pbe_=new PBEFunctional(rho.rho_);
#endif
    vxc_.resize(np_*2);
}

void PBEonGridSpin::update()
{
    get_xc_tm_.start();

    int iterative_index=rho_.getIterativeIndex();

    vector< vector<RHODTYPE> >& vrho=rho_.rho_;
//    int     ione=1;
    double  one=1.;

    Control& ct = *(Control::instance());
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();
    pb::Grid newGrid(mygrid, 2);

    pb::FDoper<RHODTYPE>* myoper_del[3];
    myoper_del[0] = new pb::Delxh4<RHODTYPE>(newGrid); 
    myoper_del[1] = new pb::Delyh4<RHODTYPE>(newGrid); 
    myoper_del[2] = new pb::Delzh4<RHODTYPE>(newGrid);
    
    pb::GridFunc<RHODTYPE>* gf_rho[2];
    for(short is=0;is<2;is++)
        gf_rho[is]=new pb::GridFunc<RHODTYPE>(&vrho[is][0],newGrid, ct.bc[0],ct.bc[1],ct.bc[2],'d');

    RHODTYPE* grad_rho[2];
    for(short is=0;is<2;is++){
        grad_rho[is]=new RHODTYPE[np_];
        memset(grad_rho[is],0,np_*sizeof(RHODTYPE));
    }

#ifdef USE_LIBXC
    double* sigma=new double[3*np_];
    memset(sigma,0,3*np_*sizeof(double));
    double* rho  =new double[2*np_];
    for ( int j = 0; j < np_; j++ ){
        for(short is=0;is<2;is++){
            rho[2*j+is]=rho_.rho_[is][j];
        }
    }
#endif


    pb::GridFunc<RHODTYPE> gf_tmp(newGrid,ct.bc[0],ct.bc[1],ct.bc[2]);
    // compute grad rho, one direction at a time
    for(short dir=0;dir<3;dir++){
        for(short is=0;is<2;is++){
            memset(grad_rho[is],0,np_*sizeof(RHODTYPE));
            //gf_rho[is]->trade_boundaries();
            myoper_del[dir]->apply(*gf_rho[is], gf_tmp);
            // convert gf_tmp back into RHODTYPE*
            gf_tmp.init_vect(grad_rho[is],'d');
        }
#ifdef USE_LIBXC
        for(short is1=0;is1<2;is1++){
            for(short is2=is1;is2<2;is2++){
                int jj=is1+is2;
                for ( int j = 0; j < np_; j++ )
                {
                    sigma[jj]+=( grad_rho[is1][j]*grad_rho[is2][j]);
                    jj+=3;
                }
            }
        }
#else
        pbe_->setGradRhoUp(dir,grad_rho[0]);
        pbe_->setGradRhoDn(dir,grad_rho[1]);
#endif
    }
    
    for(short is=0;is<2;is++)delete[] grad_rho[is];
    
    for(short i=0;i<3;i++)delete myoper_del[i];
    
    pb::GridFunc<POTDTYPE>* gf_vsigma[2];
#ifdef USE_LIBXC
    vector<double> vtmp(2*np_);
    vector<double> etmp(np_);
    vector<double> stmp(3*np_);
    xc_gga_exc_vxc(&xfunc_, np_, rho, sigma, &exc_[0], &vtmp[0], &vsigma_[0]);
    xc_gga_exc_vxc(&cfunc_, np_, rho, sigma, &etmp[0], &vxc_[0], &stmp[0]);
    delete[] sigma;
    delete[] rho;

    int nn=np_*2;
    daxpy(&nn, &one, &vxc_[0], &ione, &vtmp[0], &ione);
    nn=np_;
    daxpy(&nn, &one, &etmp[0], &ione, &exc_[0], &ione);
    
    nn=np_*3;
    daxpy(&nn, &one, &stmp[0], &ione, &vsigma_[0], &ione);
    
    //factor 2. 
    double two=2.;
    int ithree=3;
    dscal(&np_, &two, &vsigma_[0], &ithree);
    //dscal(&np_, &two, &vsigma_[np_], &ithree);
    dscal(&np_, &two, &vsigma_[2], &ithree);

    if(myspin_==0){
        for ( int j = 0; j < np_; j++ ){
            vxc_[j]=vtmp[2*j];
        }
    }else{
        for ( int j = 0; j < np_; j++ ){
            vxc_[np_+j]=vtmp[2*j+1];
        }
    }
    vector<POTDTYPE> vstmp(np_);    
    if(myspin_==0){
        for ( int j = 0; j < np_; j++ ){
            vstmp[j]=vsigma_[3*j];
        }
        gf_vsigma[0]=new pb::GridFunc<POTDTYPE>(&vstmp[0],newGrid, ct.bc[0],ct.bc[1],ct.bc[2],'d');
        for ( int j = 0; j < np_; j++ ){
            vstmp[j]=vsigma_[3*j+1];
        }
        gf_vsigma[1]=new pb::GridFunc<POTDTYPE>(&vstmp[0],newGrid, ct.bc[0],ct.bc[1],ct.bc[2],'d');
    }else{
        for ( int j = 0; j < np_; j++ ){
            vstmp[j]=vsigma_[3*j+1];
        }
        gf_vsigma[0]=new pb::GridFunc<POTDTYPE>(&vstmp[0],newGrid, ct.bc[0],ct.bc[1],ct.bc[2],'d');
        for ( int j = 0; j < np_; j++ ){
            vstmp[j]=vsigma_[3*j+2];
        }
        gf_vsigma[1]=new pb::GridFunc<POTDTYPE>(&vstmp[0],newGrid, ct.bc[0],ct.bc[1],ct.bc[2],'d');
    }
#else
    pbe_->computeXC();

    if(myspin_==0){
        assert( pbe_->pvxc1_up_!=0 );
//        Tcopy(&np_, pbe_->pvxc1_up_, &ione, &vxc_[0], &ione);
        MPcpy(&vxc_[0], pbe_->pvxc1_up_, np_);
        gf_vsigma[0]=new pb::GridFunc<POTDTYPE>(pbe_->pvxc2_upup_,newGrid, ct.bc[0],ct.bc[1],ct.bc[2],'d');
        gf_vsigma[1]=new pb::GridFunc<POTDTYPE>(pbe_->pvxc2_updn_,newGrid, ct.bc[0],ct.bc[1],ct.bc[2],'d');
    }else{
        assert( pbe_->pvxc1_dn_!=0 );
//        Tcopy(&np_, pbe_->pvxc1_dn_, &ione, &vxc_[np_], &ione);
        MPcpy(&vxc_[np_], pbe_->pvxc1_dn_, np_);
        gf_vsigma[0]=new pb::GridFunc<POTDTYPE>(pbe_->pvxc2_dnup_,newGrid, ct.bc[0],ct.bc[1],ct.bc[2],'d');
        gf_vsigma[1]=new pb::GridFunc<POTDTYPE>(pbe_->pvxc2_dndn_,newGrid, ct.bc[0],ct.bc[1],ct.bc[2],'d');
    }
    (*gf_vsigma[0])*=-1.;
    (*gf_vsigma[1])*=-1.;
#endif

    POTDTYPE* tmp=new POTDTYPE[np_];
    pb::GridFunc<POTDTYPE> gf_tmp_pot(gf_tmp);
    for(short is=0;is<2;is++){
        pb::DielFunc<POTDTYPE> diel(*gf_vsigma[is]);
        pb::PBh4<POTDTYPE> myoper(newGrid,diel); 
        // convert gf_rho to POTDTYPE
        pb::GridFunc<POTDTYPE> gf_lhs(*gf_rho[is]);
        myoper.apply(gf_lhs, gf_tmp_pot);
        // convert gf_vxc back into a POTDTYPE*
        gf_tmp_pot.init_vect(tmp,'d');
        MPaxpy(np_, one, tmp, &vxc_[np_*myspin_]);
    }

    pot_.setVxc(&vxc_[np_*myspin_],iterative_index);

    delete[] tmp;
    
    for(short isp=0;isp<2;isp++){
        delete gf_rho[isp];
        delete gf_vsigma[isp];
    }
    get_xc_tm_.stop();
}

double PBEonGridSpin::getExc()const
{    
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();

#ifdef USE_LIBXC
    double sum=0.;
    int ione=1;
    double exc = ddot(&np_, &rho_.rho_[0][0], &ione, &exc_[0], &ione);
    exc       += ddot(&np_, &rho_.rho_[1][0], &ione, &exc_[0], &ione);
    mmpi.allreduce(&exc, &sum, 1, MPI_SUM);
    exc=sum;
#else
    assert( pbe_!=NULL );
    double exc=pbe_->computeRhoDotExc();
#endif
    return exc*mygrid.vel();
}
