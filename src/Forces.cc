// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <iostream>
#include <cmath>
using namespace std;

#include "Vector3D.h"
#include "Ions.h"
#include "Potentials.h"
#include "tools.h"
#include "mdim.h"
#include "Grid.h"
#include "MGmol_blas1.h"
#include "Mesh.h"
#include "Hamiltonian.h"
#include "MPIdata.h"
#include "MGmol.h"
#include "Control.h"
#include "VariableSizeMatrix.h"
#include "DataDistribution.h"
#include "ReplicatedWorkSpace.h"
#include "KBPsiMatrixSparse.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesSparse.h"
#include "Forces.h"


#define Ry2Ha  0.5;
double  shift_R[NPTS];

Timer Forces::lforce_tm_("Forces::lforce");
Timer Forces::nlforce_tm_("Forces::nlforce");
Timer Forces::get_var_tm_("Forces::var");
Timer Forces::get_loc_proj_tm_("Forces::loc_proj");
Timer Forces::consolidate_data_("Forces::consolidate");
Timer Forces::lforce_local_tm_("Forces::lforce_local");
Timer Forces::total_tm_("Forces::total");
Timer Forces::kbpsi_tm_("Forces::KBpsi");
Timer Forces::energy_tm_("Forces::nl_energy");

double get_trilinval(const double xc, const double yc, const double zc,
                     const double h0, const double h1, const double h2,
                     const Vector3D& ref, const Vector3D& lattice,
                     RadialInter& lpot);

#if NPTS>3
double get_deriv4(double value[4]) 
{
    double sum  = (value[1] - value[0]) * 2. / (3.*DELTAC) ;
    sum -= (value[3] - value[2]) / (12.*DELTAC) ;
    return sum;
}
#endif

double get_deriv2(double value[2]) 
{
    return (value[1] - value[0]) / (2.*DELTAC) ;
}

int Forces::get_var(Ion& ion, int *pvec, double ***var_pot, 
                    double ***var_charge)
{
    get_var_tm_.start();
    
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();

    int     docount=0, incx=mygrid.dim(1) * mygrid.dim(2);
    int     ilow[3], ihi[3];

    double  rc = ion.getRC();
    assert( rc>1.e-6 );
    const double  rc2 = rc * rc;
    const double  inv_rc2=1./rc2;
    const double  rcnorm = rc*rc*rc * M_PI*sqrt(M_PI);
    const double  alpha= ion.getZion() / rcnorm;

    const int dim0 = mygrid.dim(0);
    const int dim1 = mygrid.dim(1);
    const int dim2 = mygrid.dim(2);
    
    const int gdim0 = mygrid.gdim(0);
    const int gdim1 = mygrid.gdim(1);
    const int gdim2 = mygrid.gdim(2);
    
    const double lrad=ion.getRadiusLocalPot();

    const double h0 = mygrid.hgrid(0);
    const double h1 = mygrid.hgrid(1);
    const double h2 = mygrid.hgrid(2);

    const double p0 = ion.position(0);
    const double p1 = ion.position(1);
    const double p2 = ion.position(2);

    double** potx = var_pot[0];
    double** poty = var_pot[1];
    double** potz = var_pot[2];

    double** chargex = var_charge[0];
    double** chargey = var_charge[1];
    double** chargez = var_charge[2];

    // Generate indices
    vector<vector<int> > Ai;
    Ai.resize(3);
    ion.get_Ai(Ai[0],gdim0,0);
    ion.get_Ai(Ai[1],gdim1,1);
    ion.get_Ai(Ai[2],gdim2,2);
    const int dimlx=Ai[0].size();
    const int dimly=Ai[1].size();
    const int dimlz=Ai[2].size();

    for(int i=0;i<3;i++)
    {
        ilow[i] = mygrid.istart(i);
        ihi[i]  = ilow[i] + mygrid.dim(i) - 1;
    }
    
    Vector3D ll(mygrid.ll(0),mygrid.ll(1),mygrid.ll(2));
    Vector3D position(ion.position(0),ion.position(1),ion.position(2));

    double zc = ion.lstart(2);
    const RadialInter& lpot=ion.getLocalPot();
    for(int iz = 0;iz < dimlz;iz++)
    {
        double yc = ion.lstart(1);
        for(int iy = 0;iy < dimly;iy++)
        {
            double xc = ion.lstart(0);
            for(int ix = 0;ix < dimlx;ix++)
            {
                if( (Ai[0][ix] >= ilow[0]) && (Ai[0][ix] <= ihi[0]) &&
                    (Ai[1][iy] >= ilow[1]) && (Ai[1][iy] <= ihi[1]) && 
                    (Ai[2][iz] >= ilow[2]) && (Ai[2][iz] <= ihi[2]) )
                {

                    pvec[docount] =
                         incx * (Ai[0][ix] % dim0) +
                         dim2 * (Ai[1][iy] % dim1) +
                         (Ai[2][iz] % dim2);

                    double x = xc - p0;
                    double y = yc - p1;
                    double z = zc - p2;

                    for(short ishift = 0;ishift < NPTS;ishift++)
                    {
                        double xs = x - shift_R[ishift]; // the ion moves +shift_R

                        double r2 = xs*xs + y*y + z*z;
                        double r = sqrt(r2);

                        if( r>lrad )
                        {
                            potx[ishift][docount] = 0.;
                        }
                        else
                        {
#if 0
                            potx[ishift][docount] = 
                                get_trilinval(xc-shift_R[ishift],yc,zc,
                                              h0,h1,h2,position,ll,lpot);
#else
                            potx[ishift][docount] = lpot.cubint(r);
#endif
                        }

                         chargex[ishift][docount] =
                             alpha * exp(-r2 * inv_rc2);

                    }


                    for(short ishift = 0;ishift < NPTS;ishift++)
                    {
                        double ys = y - shift_R[ishift]; // the ion moves +shift_R

                        double r2 = x*x + ys*ys + z*z;
                        double r  = sqrt(r2);

                        if( r>lrad )
                        {
                            poty[ishift][docount] = 0.;
                        }
                        else
                        {
#if 0
                            poty[ishift][docount] = 
                                get_trilinval(xc,yc-shift_R[ishift],zc,
                                              h0,h1,h2,position,ll,lpot);
#else
                            poty[ishift][docount] = lpot.cubint(r);
#endif
                        }

                        chargey[ishift][docount] =
                            alpha * exp(-r2 * inv_rc2);
                    }

                    for(short ishift = 0;ishift < NPTS;ishift++)
                    {
                         double zs = z - shift_R[ishift]; // the ion moves +shift_R

                         double r2 = x*x + y*y + zs*zs;
                         double r  = sqrt(r2);

                         if( r>lrad )
                         {
                                 potz[ishift][docount] = 0.;
                         }
                         else
                         {
#if 0
                             potz[ishift][docount] = 
                                 get_trilinval(xc,yc,zc-shift_R[ishift],
                                               h0,h1,h2,position,ll,lpot);
#else
                             potz[ishift][docount] = lpot.cubint(r);
#endif
                         }

                         chargez[ishift][docount] =
                             alpha * exp(-r2 * inv_rc2);
                    } 

                    docount++;
 
                } /* end if */

                xc += h0;

            } // end for ix
       
            yc += h1;

        } // end for iy

        zc += h2;

    } // end for iz

    assert(docount>0);
    get_var_tm_.stop();
    
    return docount;
}



void Forces::get_loc_proj(RHODTYPE *rho, const int* const pvec,
                double ***var_pot, double ***var_charge, 
                const int docount, double **loc_proj)
{
    get_loc_proj_tm_.start();
    
    Potentials& pot =hamiltonian_->potential();

    for(short dir = 0;dir < 3;dir++)
    {
        double* lproj = &(*loc_proj)[dir*NPTS];

        // pseudopotential * rho
        // - delta rhoc * vh
        for(short ishift = 0;ishift < NPTS;ishift++)
        {
            const double* const vpot =var_pot[dir][ishift];
            const double* const drhoc_ptr = var_charge[dir][ishift];
 
            for(int idx = 0;idx < docount;idx++)
            {
                const int pvidx=pvec[idx];

                lproj[ishift] += vpot[idx] * rho[pvidx];
                lproj[ishift] -= drhoc_ptr[idx]*pot.vh_rho(pvidx);  
            } 
        }
    }
    
    get_loc_proj_tm_.stop();
}

void Forces::lforce_ion(Ion& ion, RHODTYPE *rho, double **loc_proj)
{
    double  ***var_pot, ***var_charge;

    Mesh* mymesh = Mesh::instance();
    const int numpt=mymesh->numpt();    

    DIM3(var_pot, 3, NPTS, numpt, double);
    int* pvec=new int[numpt];
    DIM3(var_charge, 3, NPTS, numpt, double);

    // generate pvec, var_pot and var_charge for this ion 
    int docount=get_var(ion, pvec, var_pot, var_charge);
    
    get_loc_proj(rho, pvec, var_pot, var_charge, docount, loc_proj);


    D3FREE(var_charge);
    delete[] pvec;
    D3FREE(var_pot);
}




void Forces::lforce(Ions& ions, RHODTYPE *rho)
{
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();
//    Control& ct = *(Control::instance());    

    lforce_tm_.start();
    
//    if(ct.short_sighted)
//    {
       //cout<<"max Vl radius = "<<ions.getMaxVlRadius()<<endl;
       
       int buffer_size = 3*NPTS;
       double* loc_proj = new double[buffer_size];
       
       lforce_local_tm_.start();
       
       int cols[buffer_size];
       for(int i=0; i<buffer_size; i++)
          cols[i] = i;
       
       VariableSizeMatrix<sparserow> loc_proj_mat("locProj",ions.overlappingVL_ions().size());
       // Loop over ions with potential overlaping with local subdomain
       vector<Ion*>::const_iterator ion=ions.overlappingVL_ions().begin();
       while(ion!=ions.overlappingVL_ions().end())
       {
           int index=(*ion)->index();
           for(short dir=0;dir<3*NPTS;dir++)loc_proj[dir]=0.;
           lforce_ion(**ion, rho, &loc_proj);

           /* insert row into 2D matrix */
           loc_proj_mat.insertNewRow(buffer_size, index, &cols[0], loc_proj, true);
           
           ion++;
       } 
       
       lforce_local_tm_.stop();
       consolidate_data_.start();
       
       /* Distribute/ gather data */
       /* ions have rectangular domain and we only care about gathering on 
        * ions centered on the local domain, hence we set append=false
       */
       const pb::PEenv& myPEenv=mymesh->peenv();
       double domain[3]={mygrid.ll(0),mygrid.ll(1),mygrid.ll(2)};
       DataDistribution distributor("lforce",ions.getMaxVlRadius(), myPEenv, domain );
       //the first time through, we may need to recompute the buffer size for
       //sparse matrix communications since we are dealing with a matrix not seen before
       static bool first_time=true;
       if( first_time )
       {
           DataDistribution::enforceComputeMaxDataSize();
           first_time=false;
       }
       distributor.augmentLocalData(loc_proj_mat, false); 
       
       consolidate_data_.stop();

       vector<Ion*>::const_iterator lion=ions.local_ions().begin();
       while(lion!=ions.local_ions().end())
       {
           // Forces opposed to the gradient
           int index=(*lion)->index();
           int *rindex = (int *)loc_proj_mat.getTableValue(index);
           assert(rindex != NULL);
           memset(loc_proj, 0, buffer_size*sizeof(double));
           loc_proj_mat.row_daxpy(*rindex, buffer_size, mygrid.vel(), loc_proj);

           (*lion)->add_force(-get_deriv2(loc_proj),
                              -get_deriv2((loc_proj + NPTS)),
                              -get_deriv2((loc_proj + 2*NPTS)) );
           
           lion++;
       }
       delete [] loc_proj;
//    }    
//    else
/*
    {
       double  **loc_proj;
    
       DIM2(loc_proj, ions.getNumIons(), 3*NPTS, double);

       init_loc_proj(loc_proj,ions.getNumIons());

       // Loop over ions 
       vector<Ion*>::const_iterator ion=ions.overlappingVL_ions().begin();
       while(ion!=ions.overlappingVL_ions().end()){
           
           int index=(*ion)->index();
           lforce_ion(**ion, rho, &loc_proj[index]);
           
           ion++;
       } 

       int      n = 3*NPTS*ions.getNumIons();
       pb::my_dscal(n,mygrid.vel(),*loc_proj);
       global_sums_double(*loc_proj, n);

       vector<Ion*>::iterator lion=ions.local_ions().begin();
       while(lion!=ions.local_ions().end()){
           // Forces opposed to the gradient
           int index=(*lion)->index();
           (*lion)->add_force(-get_deriv2(&loc_proj[index][0]),
                              -get_deriv2(&loc_proj[index][1*NPTS]),
                              -get_deriv2(&loc_proj[index][2*NPTS]) );
           
           lion++;
       }
       D2FREE(loc_proj);
    }
*/
#ifdef HAVE_TRICUBIC
    Potentials& pot =hamiltonian_->potential();
    if( pot.withVext() )
    {
        double position[3];
        double grad[3];
        vector<Ion*>::iterator ion=ions.local_ions().begin();
        while(ion!=ions.local_ions().end())
        {
            (*ion)->getPosition(position);

            pot.getGradVext(position,grad);

            const double charge=(*ion)->getZion();

            (*ion)->add_force(grad[0]*charge,
                              grad[1]*charge, 
                              grad[2]*charge);
         //if( onpe0 )     
         //(*MPIdata::sout)<<"External force on Ion "<<(*ion)->name()<<": "
         //                              <<grad[0]*charge<<","
         //                              <<grad[1]*charge<<","
         //                              <<grad[2]*charge<<endl;
            ion++; 
        }
    }
#endif
    
    lforce_tm_.stop();
}

// Get the nl energy as the trace of loc_kbpsi*mat_X for several loc_kbpsi
// result added to erg

void Forces::nlforceSparse(LocGridOrbitals& orbitals, Ions& ions)
{
    if( ions.getNumIons()==0 )
        return;

    Control& ct = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    // first check if any NL forces computation necessary
    if( !ions.hasNLprojectors() )
    {
        if( onpe0 )cout<<"No nl forces!!"<<endl;
        return;
    }

    nlforce_tm_.start();
    
    kbpsi_tm_.start();
    KBPsiMatrixSparse*** kbpsi=new KBPsiMatrixSparse**[3];
    
    // compute all kbpsi matrices for all shifts
    for(short dir=0;dir<3;dir++)
    {
        kbpsi[dir] = new KBPsiMatrixSparse*[NPTS];
        for(int npt=0;npt<NPTS;npt++)
        {
            kbpsi[dir][npt]=new KBPsiMatrixSparse(NULL,false);

            double shift[3]={0.,0.,0.};
            shift[dir]=shift_R[npt];
            Ions shifted_ions(ions, shift); ///***
            
            kbpsi[dir][npt]->setup(shifted_ions, orbitals);
            kbpsi[dir][npt]->computeAll(shifted_ions, orbitals);
        }
    }
    kbpsi_tm_.stop();
    
    energy_tm_.start();
    map<int,double*> erg;
    if( ct.short_sighted )
    {
        ProjectedMatricesSparse* projmatrices =
            dynamic_cast<ProjectedMatricesSparse*>(proj_matrices_);
        assert( projmatrices );
        DensityMatrixSparse& dm( projmatrices->getDM());

        // loop over all the ions
        // parallelization over ions by including only those centered in subdomain
        const vector<Ion*>::const_iterator iend=ions.local_ions().end();
        vector<Ion*>::const_iterator       ion =ions.local_ions().begin();
        while(ion!=iend)
        {
            vector<int> gids;
            (*ion)->getGidsNLprojs(gids);
            vector<short> kbsigns;
            (*ion)->getKBsigns(kbsigns);
            
            const short nprojs=(short)gids.size();
            for(short i=0;i<nprojs;i++)
            {
                const int gid=gids[i];
                
                double* zeros=new double[3*NPTS];
                memset(zeros,0,3*NPTS*sizeof(double));
                erg.insert(pair<int,double*>(gid,zeros));
            }
            
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(short ii=0;ii<nprojs*3*NPTS;ii++)
            {
                short ip=ii/(3*NPTS);
                const int gid=gids[ip];
                const double kbmult=(double)kbsigns[ip];
                
                short it = ii%(3*NPTS);
                //for(short it=0;it<3*NPTS;it++)
                {
                    short dir   =it/NPTS;
                    short ishift=it%NPTS;
                    double alpha=kbpsi[dir][ishift]->getTraceDM(gid, dm);
                    erg[gid][NPTS*dir+ishift]=alpha*kbmult;
                }
            }
            ion++;
        }      
    }else{
        ReplicatedWorkSpace& wspace( ReplicatedWorkSpace::instance() );
        const int ndim = wspace.getDim();
        DISTMATDTYPE* work_DM_matrix=wspace.square_matrix();

        assert( dynamic_cast<ProjectedMatrices*>( proj_matrices_ ) );
        ProjectedMatrices* projmatrices=dynamic_cast<ProjectedMatrices*>( proj_matrices_ );
        projmatrices->getReplicatedDM(work_DM_matrix);

        // loop over all the ions
        // parallelization over ions by including only those centered in subdomain
        const vector<Ion*>::const_iterator iend=ions.local_ions().end();
        vector<Ion*>::const_iterator       ion =ions.local_ions().begin();
        while(ion!=iend)
        {
            vector<int> gids;
            (*ion)->getGidsNLprojs(gids);
            vector<short> kbsigns;
            (*ion)->getKBsigns(kbsigns);
            
            const short nprojs=(short)gids.size();
            for(short i=0;i<nprojs;i++)
            {
                const int gid=gids[i];
                const double kbmult=(double)kbsigns[i];
                
                double* zeros=new double[3*NPTS];
                memset(zeros,0,3*NPTS*sizeof(double));
                erg.insert(pair<int,double*>(gid,zeros));
                
                for(short dir=0;dir<3;dir++)
                for(short ishift = 0;ishift < NPTS;ishift++)
                {
                    double alpha=kbpsi[dir][ishift]->getTraceDM(gid, work_DM_matrix, ndim);
                    erg[gid][NPTS*dir+ishift]=alpha*kbmult;
                }
            }
            ion++;
        }
    }
    energy_tm_.stop();
    
    
    // release memory
    for(short dir=0;dir<3;dir++)
    {
        for(short npt=0;npt<NPTS;npt++)
        {
            delete kbpsi[dir][npt];
        }
        delete[] kbpsi[dir];
    }
    delete[] kbpsi;
    
    // compute forces on each ion by finite differences
    const double factor=-1.*Ry2Ha;
    vector<Ion*>::iterator             iion = ions.local_ions().begin();
    const vector<Ion*>::const_iterator iend = ions.local_ions().end();
    while(iion!=iend)
    {
        vector<int> gids;
        (*iion)->getGidsNLprojs(gids);
        
        const short nprojs=(short)gids.size();
        for(short i=0;i<nprojs;i++)
        {
            const int gid=gids[i];
            
            
            double  ff[3]={get_deriv2(&erg[gid][NPTS*0])*factor,
                           get_deriv2(&erg[gid][NPTS*1])*factor,
                           get_deriv2(&erg[gid][NPTS*2])*factor};
            
            if(mmpi.nspin()==2)
            {
                double sum[3]={0.,0.,0.};
                mmpi.allreduceSpin(&ff[0],&sum[0],3,MPI_SUM);
                for(short dir=0;dir<3;dir++)ff[dir]=sum[dir];
            }

            (*iion)->add_force(ff[0], ff[1], ff[2]);
        }
        
        ++iion;
    }
    
    map<int,double*>::iterator ierg=erg.begin();
    while( ierg!=erg.end() )
    {
        delete[] ierg->second;
        ++ierg;
    }

    nlforce_tm_.stop();
}

void Forces::force(LocGridOrbitals& orbitals, Ions& ions)
{
#ifdef USE_BARRIERS
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.barrier();       
#endif
    
    total_tm_.start();
    
    const int numpt=rho_->rho_[0].size();
    double one=1.;
    
    vector<RHODTYPE> rho_tmp;
    if(rho_->rho_.size()>1)
    {
        rho_tmp.resize(numpt);
        
        memcpy(&rho_tmp[0],&(rho_->rho_[0][0]),numpt*sizeof(RHODTYPE));
        MPaxpy(numpt,one,&rho_->rho_[1][0],&rho_tmp[0]);
    }
    
    vector<RHODTYPE>& rho = (rho_->rho_.size()>1) ? rho_tmp : rho_->rho_[0];
    
    shift_R[0] = -DELTAC;
    shift_R[1] =  DELTAC;
#if NPTS>3
    shift_R[2] = -2.*DELTAC;
    shift_R[3] =  2.*DELTAC;
#endif
    Control& ct = *(Control::instance());

    // Zero out forces 
    ions.resetForces();

    // Get the ion-ion component and store.
    ions.iiforce( ct.bcPoisson );

    // Add the non-local forces
#ifdef USE_BARRIERS
    mmpi.barrier();       
#endif
//    if(ct.short_sighted)
        nlforceSparse(orbitals,ions);
//    else
//        nlforce(orbitals,ions);

    // Add the local forces
#ifdef USE_BARRIERS
    mmpi.barrier();       
#endif
    lforce(ions,&rho[0]);
    
#ifdef USE_BARRIERS
    mmpi.barrier();       
#endif
    
    total_tm_.stop();
}

