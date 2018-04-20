// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$

//Old interface functions.

#include <iostream>
#include <string>
#include <cassert>
#include <list>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "hdf5.h"

#include "Control.h"
#include "Mesh.h"

#include "Hamiltonian.h"
#include "Species.h"
#include "Ions.h"
#include "HDFrestart.h"
#include "Potentials.h"

#include "GridFactory.h"
#include "KBPsiMatrix.h"
#include "ConstraintSet.h"
#include "LocalizationRegions.h"
#include "tools.h"
#include "BlockVector.h"
#include "MGmol.h"

#define max(a,b) (((a)<(b)) ? (b) : (a))

//#define DEBUG 1

// Description of the run 
extern string description;

//vector<double> occupation;

const double ry2ev=13.605804;


int MGmol::read_params1(ifstream* tfile)
{
    Control& ct = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    
    if( mmpi.instancePE0() )
    {
        string query;

        // Read Laplacian operator type
        (*tfile)>>ct.lap_type;
        read_comments(*tfile);

        // Read in the boundary condition flag for Vh
        string queryBcPoisson;
        if ( !getline ( *tfile, queryBcPoisson ) ) {
            (*MPIdata::serr)<<"readParameters: cannot read BC"<<endl;
            return -1;
        }
        stringstream ssBcPoisson ( queryBcPoisson );
        for(int i=0;i<3;i++)
            ssBcPoisson>>ct.bcPoisson[i];            
        if( ct.bcPoisson[0] == 2 
         || ct.bcPoisson[1] == 2
         || ct.bcPoisson[2] == 2 ){
            if( !(ssBcPoisson>>ct.multipole_order) )
                ct.multipole_order=1; // default value
        }

        for(int i=0;i<3;i++)
        if ((ct.bcPoisson[i] != 0)
         && (ct.bcPoisson[i] != 1)
         && (ct.bcPoisson[i] != 2) ){
            cerr<<"readParameters: invalid boundary conditions"<<endl;
            return -1;
        }
        if(onpe0)os_<<"Poisson BC: "<<ct.bcPoisson[0]<<","
                                    <<ct.bcPoisson[1]<<","
                                    <<ct.bcPoisson[2]<<endl;

        read_comments(*tfile);
 
        // Dielectric or not?
        short diel;
        (*tfile)>>diel;
        if( (diel%10)!=0 && (diel%10)!=1 ){
            cerr<<"Flag diel should be 0, 1, 10, or 11"<<endl;
            return -1;
        }   
        // # of Hartree iterations
        (*tfile)>>ct.vh_init;
        if( ct.vh_init<0 || ct.vh_init>100 ){
            cerr<<"Invalid parameter ct.vh_init"<<endl;
            return -1;
        }
        (*tfile)>>ct.vh_its;
        if( ct.vh_its<0 || ct.vh_its>100 ){
            cerr<<"Invalid parameter ct.vh_its"<<endl;
            return -1;
        }
        // MG preconditioner parameters for PCG Poisson solver
        short nu1=2;
        (*tfile)>>nu1;
        if(tfile->fail())
        {
            tfile->clear();
        }        
        short nu2=2;
        (*tfile)>>nu2;
        if(tfile->fail())
        {
            tfile->clear();
        } 
        short nlev=10;
        (*tfile)>>nlev;
        if(tfile->fail())
        {
            tfile->clear();
        } 
        ct.setDielAlgo(diel, nu1, nu2, nlev);
        
        read_comments(*tfile);

        // dielectric parameters
        (*tfile)>>ct.rho0;
        if(ct.rho0>.1){
            cerr<<"rho0="<<ct.rho0<<" is too large"<<endl;
            return -1;
        }

        (*tfile)>>ct.drho0;
        read_comments(*tfile);
        if(onpe0)os_<<"ct.rho0="<<ct.rho0<<", ct.drho0="<<ct.drho0<<endl;

        //bias
        (*tfile)>>ct.screening_const;
        read_comments(*tfile);
 
        // Short-sighted algorithm 
        char short_sighted;
        float spread_factor, fgmres_tol, ilu_droptol;
        short fgmres_kim, fgmres_maxits, ilu_lof, ilu_maxfil, ilu_type;
        (*tfile)>>short_sighted;
        if( short_sighted=='p' || short_sighted=='d' ){
            if(onpe0)os_<<"WARNING: Mix type is an obsolete option."<<endl;
            ct.short_sighted = 0;
        }else{
            if( short_sighted=='s' ){
                ct.short_sighted = 1;
                if(onpe0)os_<<"Short-sighted algorithm selected..."<<endl;
                
                /* solver parameters for computing inverse entries */ 
                (*tfile)>>spread_factor>>fgmres_tol>>fgmres_kim>>fgmres_maxits>>ilu_type>>ilu_droptol>>ilu_lof>>ilu_maxfil;
                int info=ct.setShortSightedSolverParameters(spread_factor, fgmres_tol, ilu_droptol, fgmres_kim, fgmres_maxits, ilu_lof, ilu_maxfil, ilu_type);
                if( info==-1)return -1;
            }
            if( short_sighted=='l' )
                ct.short_sighted = 0;
        }
        if( ct.short_sighted!=0 && ct.short_sighted!=1 ){
            cerr<<"invalid parameter for short-sighted option: "<<ct.short_sighted<<endl;
            return -1;
        }
    
        read_comments(*tfile);

        // Algorithm type for wavefunctions optimization
        (*tfile)>>ct.it_algo_type;
        if( ct.it_algo_type==2 )(*tfile)>>ct.dm_inner_steps;
        short coloring_algo=0; // 0=global RLF
        (*tfile)>>coloring_algo;
        if(tfile->fail())
        {
            tfile->clear();
            coloring_algo=0;
        }
        ct.setColoringAlgo(coloring_algo);
        read_comments(*tfile);

        // Mixing parameters
        (*tfile)>>ct.mix_pot>>ct.dm_mix;
        if( ct.mix_pot>2. || ct.dm_mix>2. || ct.mix_pot<0. || ct.dm_mix<0.){
            cerr<<"Invalid mixing parameters"<<endl;
            return -1;
        }
        if( ct.it_algo_type>1 )ct.mix_pot=1.;
        read_comments(*tfile);

        (*tfile)>>ct.atoms_dyn;
        if( ct.atoms_dyn!=0 && ct.atoms_dyn!=2 && ct.atoms_dyn!=6 && ct.atoms_dyn!=7 ){
            cerr<<ct.atoms_dyn<<" is invalid parameter for md method"<<endl;
            return -1;
        }
        if(onpe0)os_<<"parameter for md method: "<<ct.atoms_dyn<<endl;

        //unused option
        char dummy;
        (*tfile)>>dummy;
        read_comments(*tfile);

        {
            // Number of steps 
            if ( !getline ( *tfile, query ) ) {
                cerr<<"Cannot read number of step!!!"<<endl;
                return -1;
            }
            stringstream ss ( query );
            ss>>ct.num_MD_steps;
            os_<<"Run "<<ct.num_MD_steps<<" iterations"<<endl;
            if( ct.atoms_dyn==2 || ct.atoms_dyn==6 || ct.atoms_dyn==7 ){
                ss>>ct.max_electronic_steps;
                if( ct.max_electronic_steps<1 ){
                    cerr<<"Run method requires a positive number of inner iterations"
                        <<endl;
                    return -1;
                }
                os_<<"Run "<<ct.max_electronic_steps<<" inner iterations"<<endl;
                
                if( !(ss>>ct.md_print_freq) )
                    ct.md_print_freq=1;
                if( !(ss>>ct.md_print_filename) )
                    ct.md_print_filename="0";
            }else{
                ct.max_electronic_steps=ct.num_MD_steps;
            }
        }
        read_comments(*tfile);
 
        // Localization centers updates
        int info=ct.readLRupdateInfo(tfile);
        if( info<0 )return -1;
        read_comments(*tfile);

        // Type of MD: thermostat? Temperature? Relaxation time?
        ct.readThermostatInfo(tfile);
        read_comments(*tfile);

        //(*tfile)>>ct.pdamp; // not used anymore
        //ct.pdamp=0.;
        //use kernel functions for projections for AOMM approach?
        if ( !getline ( *tfile, query ) ) {
            cerr<<"Cannot read use_kernel_functions flag!!!"<<endl;
            return -1;
        }
        stringstream ssk ( query );
        if( !(ssk>>ct.use_kernel_functions) )
            ct.use_kernel_functions=0; // default value

        if( ct.use_kernel_functions )
        {
            cout<<"AOMM algorithm..."<<endl;
        }
        read_comments(*tfile);

        // Quench method 
        (*tfile)>>ct.wf_dyn;
        if( ct.wf_dyn!=0 && ct.wf_dyn!=1 ){
            cerr<<ct.wf_dyn<<" is an invalid quench option!!!"<<endl;
            return -1;
        }
        if( ct.wf_dyn==1 ){
           (*tfile)>>ct.wf_m;
           (*tfile)>>ct.betaAnderson;
            if(onpe0)os_<<"quench: Anderson extrapolation with m="<<ct.wf_m
                <<" and beta="<<ct.betaAnderson<<endl;
        }
        read_comments(*tfile);
 
        (*tfile)>>ct.wf_extrapolation;
        if( !(ct.wf_extrapolation==0
           || ct.wf_extrapolation==1
           || ct.wf_extrapolation==2) ){
            cerr<<"Invalid option for WF extrapolation in MD"<<endl;
            return -1;
        }
        if( ct.wf_extrapolation==0 )
        {
            if(onpe0)os_<<"WF extrapolation: reversible scheme"<<endl;
        }
        if( ct.wf_extrapolation==1 )
            if(onpe0)os_<<"WF extrapolation: 2nd order scheme"<<endl;
        if( ct.wf_extrapolation==2 )
            if(onpe0)os_<<"WF extrapolation: 3rd order scheme"<<endl;
        read_comments(*tfile);

        (*tfile)>>ct.iprint_residual;
        //if( ct.iprint_residual>0 )
        //    os_<<"Print residuals every "<<ct.iprint_residual
        //        <<" iterations"<<endl;
        read_comments(*tfile);

        // Ionic timestep
        if ( !getline ( *tfile, query ) ) {
            cerr<<"Cannot read time step!!!"<<endl;
            return -1;
        }
        stringstream ss ( query );
        ss>>ct.dt;            
        if( ct.atoms_dyn==2 && ct.dt>0. ){
            if( !(ss>>ct.enforceVmass0) )
                ct.enforceVmass0=0; // default value
        }
        //os_<<"ct.dt="<<ct.dt<<endl;
        //os_<<"ct.enforceVmass0="<<ct.enforceVmass0<<endl;
        read_comments(*tfile);
 
        // convergence criterion 
        if ( !getline ( *tfile, query ) ) {
            cerr<<"Cannot read convergence criterion!!!"<<endl;
            return -1;
        }
        stringstream sstol ( query );
        sstol>>ct.conv_tol;
        sstol>>ct.tol_forces;  
        if( !(sstol>>ct.conv_tol_stop) )
            ct.conv_tol_stop=1.; // default value
        //(*tfile)>>ct.conv_tol>>ct.tol_forces;
        read_comments(*tfile);
 
        // Exchange correlation potential type flag 
        (*tfile)>>ct.xctype;
        switch(ct.xctype){
            case 0:
#ifdef DEBUG
                if(onpe0)os_<<"LDA"<<endl;
#endif
                break;
            case 2:
#ifdef DEBUG
                if(onpe0)os_<<"PBE"<<endl;
#endif
                break;
            default:
                cerr<<"EX and CORR undefined!!!"<<endl;
                return -1;
        }
        read_comments(*tfile);

        // Checkpoint count 
        (*tfile)>>ct.checkpoint;
        read_comments(*tfile);
   
        (*tfile)>>ct.verbose;
#ifdef PRINT_OPERATIONS
        ct.verbose=max(ct.verbose,2);
#endif    
        read_comments(*tfile);
   
        (*tfile)>>ct.occ_width;    
#ifdef DEBUG
        if(onpe0)os_<<"Occupation width = "<<ct.occ_width<<"[eV]"<<endl;
#endif
        ct.occ_width /= ry2ev; // convert to Ry
        read_comments(*tfile);
 
        // Type of orbitals
        (*tfile)>>ct.orbital_type;
        (*tfile)>>ct.dot_product_type;
        if(tfile->fail())
        {
            ct.dot_product_type=0; // default value
            tfile->clear();
        }

        switch(ct.orbital_type)
        {
            case 0:
                if(onpe0)os_<<"Eigenfunctions"<<endl;
                break;
            case 1:
                if(onpe0)os_<<"Nonorthogonal orbitals with dot product type "<<ct.dot_product_type<<endl;
                break;
            case 2:
                if(onpe0)os_<<"Orthonormal orbitals with dot product type "<<ct.dot_product_type<<endl;
                break;
            default:
                cerr<<"Orbitals type undefined!!!"<<endl;
                return -1;
        }
        read_comments(*tfile);
 
        // Initialization with localized orbitals (1) or not (0)
        (*tfile)>>ct.init_loc>>ct.init_type;
        if( ct.init_type==1 )
        {
            (*tfile)>>ct.init_rc;
            if(onpe0)os_<<"Initialize orbitals with Gaussians of width "<<ct.init_rc<<endl;
        }
        else if( ct.init_type==2 )
        {
            if(onpe0)os_<<"Initialize orbitals with Fourier basis"<<endl;
        }
        else if( ct.init_type==0 )
        {
            if(onpe0)os_<<"Initialize orbitals with random values"<<endl;
            if( ct.init_loc )
            {
                (*tfile)>>ct.init_rc;
                if(onpe0)os_<<"Initialize orbitals of radius "<<ct.init_rc<<endl;
            }
        }
        else
        {
            cerr<<"Invalid initial type for  orbitals!!!"<<endl;
            return -1;
        }

        read_comments(*tfile);
 
        // MD Integration flag 
        int mdflag;
        (*tfile)>>mdflag;
        if( mdflag!=0 ){
            cerr<<mdflag<<" is an invalid time integration method!!!"<<endl;
            return -1;
        }
        read_comments(*tfile);
 
        // Number of states 
        (*tfile)>>ct.numst;
        if( ct.numst<0 ){
            cerr<<"Invalid number of states"<<endl;
            return -1;
        }
        read_comments(*tfile);

#ifdef DEBUG
        os_<<ct.numst<<" states\n";
#endif        
    }
    return 0;
}

//old interface
int MGmol::readParameters(ifstream* tfile, bool& cell_relative)
{ 
    Control& ct = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    
    // cell information
    double origin[3];
    float end[3];
 
    if( mmpi.instancePE0() )
    {
        // Read in the description 
        read_comments(*tfile);
        (*tfile)>>description;
        read_comments(*tfile);
        
        // Lattice constants 
        (*tfile)>>origin[0]>>origin[1]>>origin[2];
        (*tfile)>>end[0]>>end[1]>>end[2];
        read_comments(*tfile);
 
        // Grid size 
        (*tfile)>>ct.ngpts_[0]>>ct.ngpts_[1]>>ct.ngpts_[2];
        read_comments(*tfile);
 
    }
    mmpi.bcast(&origin[0], 3);
    mmpi.bcast(&end[0],    3);
    int ngpts[3]={(int)ct.ngpts_[0],(int)ct.ngpts_[1],(int)ct.ngpts_[2]};
    mmpi.bcast(&ngpts[0],  3);

    const double cell[3]={end[0]-origin[0],end[1]-origin[1],end[2]-origin[2]};

    if(onpe0){
        os_<<" cell=("<<cell[0]<<","<<cell[1]<<","<<cell[2]<<")"<<endl;
        os_<<" grid=("<<ct.ngpts_[0]<<","<<ct.ngpts_[1]<<","<<ct.ngpts_[2]<<")"<<endl;
    }

    if( ct.ngpts_[0]<1 )return -1;
    if( ct.ngpts_[1]<1 )return -1;
    if( ct.ngpts_[2]<1 )return -1;

    if( mmpi.instancePE0() ){
        // Number of species 
        (*tfile)>>ct.num_species;
        read_comments(*tfile);
        if(onpe0)os_<<ct.num_species<<" species"<<endl;
    }
#ifdef USE_MPI
    mmpi.bcast(&ct.num_species);
#endif    
    if( ct.num_species<0 )return -1;
    if( ct.num_species==0 ){
        if(onpe0)os_<<"Warning: no atoms in cell"<<endl;
    }
    
    ct.readPotFilenames(tfile);

    ct.readRestartInfo(tfile);
    
    if( mmpi.instancePE0() )read_comments(*tfile);

    ct.readRestartOutputInfo(tfile);
    
    if( mmpi.instancePE0() )read_comments(*tfile);
 
    int flag=read_params1(tfile);

    ct.sync();
    
    Mesh::setup(comm_,ct.ngpts_,origin,cell,ct.lap_type);
       
    hamiltonian_ = new Hamiltonian();
    Potentials& pot =hamiltonian_->potential();

    ct.registerPotentials(pot);
/*
    if( ct.wf_extrapolation==0 )
    {
        BlockVector<ORBDTYPE>::incMaxAllocInstances(3);
    }
*/

    mmpi.bcast(&flag);

    if( flag<0 )return flag;

    Mesh* mymesh = Mesh::instance();
    const pb::Grid&  myGrid  = mymesh->grid();

    assert( ct.xctype==0 ||  ct.xctype==2 );
    assert( ct.mix_pot < 2. && ct.mix_pot>0. );

    // Occupations of the states 
    ct.readOccupations(tfile);
    
    int cflag;
    int dummy;
    if( mmpi.instancePE0() )
    {
         // Number of ions 
        (*tfile)>>dummy; // ignored, obsolete

        read_comments(*tfile);

        int precond_type, mg_levels, project_out_psd;
        float precond_factor;
 
        (*tfile)>>precond_type;
        
        // Wavefunction smoothing timestep for the subiteration 
        (*tfile)>>precond_factor>>mg_levels>>project_out_psd;
        int info=ct.setPreconditionerParameters(precond_type,precond_factor,
                                                (project_out_psd==1),mg_levels,
                                                myGrid.hmax());
        if( info==-1)return -1;
        read_comments(*tfile);
 
        // Number of iteration without potential change in a quench 
        (*tfile)>>ct.max_changes_pot;
        read_comments(*tfile);
 
        // Cutting radius 
        float cut_radius;
        (*tfile)>>cut_radius>>ct.orthof;
        
        float tol_eigen_gram=0.;
        float min_distance_centers=0.;
        (*tfile)>>tol_eigen_gram;
        if(tfile->fail())
        {
            tfile->clear();
            tol_eigen_gram=0.;
            if(onpe0)os_<<"uses default value tol_eigen_gram = "<<tol_eigen_gram<<endl;
        }else{
        
            (*tfile)>>min_distance_centers;
            if(tfile->fail())
            {
                tfile->clear();
                min_distance_centers=0.;
            }
        }
        ct.setLocMode(cut_radius, myGrid.ll(0),myGrid.ll(1),myGrid.ll(2),
                      min_distance_centers);
        ct.setTolEigenvalueGram(tol_eigen_gram);

#ifdef DEBUG
        os_<<" ct.orthof="<<ct.orthof<<endl;
#endif        
        read_comments(*tfile);


         // Absolute or cell relative coordinates 
        // 0=cell relative and 1=absolute 
        // override type, positions and force control characters 
         // 0=do not override and 1=override 
        (*tfile)>>cflag>>ct.override_restart;
        if( !( cflag==0 || cflag==1 ) ){
            cerr<<"Invalid coordinate system"<<endl;
            return -1;
        }
 
        read_comments(*tfile);
 
         // Read in flag to compute Wannier centers
        (*tfile)>>ct.wannier_transform_type;
        assert( ct.wannier_transform_type==0 
             || ct.wannier_transform_type==1 
             || ct.wannier_transform_type==2 );
        //if( ct.lr_update )ct.wannier_transform_type=1;
        read_comments(*tfile);
 
        // transfer matrices flag 
        (*tfile)>>ct.tmatrices;
        assert( ct.tmatrices==1 || ct.tmatrices==0 );
        
        if(ct.override_restart)os_<<"override coordinates"<<endl;

        read_comments(*tfile); 
    }

#ifdef USE_MPI
#ifdef DEBUG
    if(onpe0)
        os_<<"MPI Bcast in readParameters!!!"<<endl;
#endif
    mmpi.bcast(&cflag);
#endif        
    ct.sync();
    cell_relative = !(bool)cflag;

    mymesh->subdivGridx( ct.getMGlevels() );

#ifdef DEBUG
    if(onpe0)
        os_<<"subdivx="<<mymesh->subdivx()<<endl;
#endif
    assert(cell[0]>0.);
    assert(cell[1]>0.);
    assert(cell[2]>0.);
     

#ifdef DEBUG
    if(onpe0)
        os_<<" Create "<<ct.num_species<<" species"<<endl;
#endif

    ct.setSpecies(pot);
    
#ifdef DEBUG
    printWithTimeStamp(" Pseudopotential initialized",os_);
#endif
    mmpi.barrier();
    
    ct.checkNLrange();
    
    if( ct.diel )
        pot.turnOnDiel();

    if( ct.screening_const>0. && ct.lap_type != 0){
        ct.screening_const=0.;
        if(onpe0)
            os_<<"Coulomb screening not used (not implemented for this FD scheme)"
                <<endl;
    }
    
    //set parameter values not available in that interface
    ct.setDefaultValues();
    
    ct.adjust();
    
    return ct.checkState();
}

int MGmol::readLRsFromInput(ifstream* tfile)
{ 
    assert( lrs_!=0 );

    Control& ct ( *(Control::instance()) );
    
    if( ct.verbose>0 )
        printWithTimeStamp("readLRsFromInput",os_);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    if( ct.restart_info>2 )return 0;

    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();

    const double lattice[3]={mygrid.ll(0),mygrid.ll(1),mygrid.ll(2)};
    const double origin[3] ={mygrid.origin(0), mygrid.origin(1), mygrid.origin(2)};
    const double end[3] ={mygrid.origin(0)+mygrid.ll(0), 
                          mygrid.origin(1)+mygrid.ll(1), 
                          mygrid.origin(2)+mygrid.ll(2)};
    
    const bool read_radius = ct.adaptiveLRsizes() ? true : false;

    if( mmpi.instancePE0() )
    { 
        if( ct.verbose>0 )os_<<"Trying to read "<<ct.numst <<" orbital centers from input file..."<<endl;
        const vector<Ion*>& list_ions(ions_->list_ions());
        vector<Ion*>::const_iterator ion=list_ions.begin();
        
        if( tfile!=0 )read_comments(*tfile);
        int count=0;
        
        // read as many centers as there are functions
        for( int i=0;i<ct.numst;i++ )
        {
            double crds[3]={0.,0.,0.};
            bool flag=false;
            
            if( tfile!=0 )
            for(short j=0;j<3;j++)
            {
                string sread;
                (*tfile)>>sread;
                if( tfile->fail() )
                {
                    os_<<"WARNING: Failed reading Localization center... Use atomic position for center "<<i<<endl;
                    flag=false;
                    break;
                }
                else
                {
                    flag=true;
                    crds[j]=atof(sread.c_str());
                }
            }
            
            // set center to ionic position if not read from input file
            if( !flag )
            {
                for(int k=0;k<3;k++)crds[k]=(*ion)->position(k);
                ion++;
                if( ion==list_ions.end() )ion=list_ions.begin();
            }
            
            // move coordinates inside domain boundaries
            for(int j=0;j<3;j++)
            {
                while( crds[j]<origin[j] ) crds[j]+=lattice[j];
                while( crds[j]>=end[j]   ) crds[j]-=lattice[j];
                
                assert(crds[j]>-1000.);
                assert(crds[j]< 1000.);
            }
            
            float radius=ct.cut_radius;
            if( flag && ct.isLocMode() && read_radius )
            {
                float tmp;
                (*tfile)>>tmp;
                if( tfile->fail() )
                {
                    os_<<"WARNING: Failed reading Localization radius!!! Use uniform radius..."<<endl;
                }
                else
                {
                    radius=tmp;
                }
            }
            
            //this is where a new function with a new gid is created
            Vector3D tmp_center(crds[0], crds[1], crds[2]);
            lrs_->push_back_global(tmp_center,radius);
            if ( ct.verbose>2 )
                (*MPIdata::sout)<<"Added LR with center "<<tmp_center<<" and radius "<<radius<<endl;
            count++;
            
            // finish reading line
            if( flag )while( tfile->get()!='\n');
#ifdef DEBUG
            os_<<" Read orbital center ("
               <<crds[0]<<","
               <<crds[1]<<","
               <<crds[2]<<")"<<endl;
#endif 
        }
        
        if ( ct.verbose>1 && onpe0 ) os_<<"Randomize gids..."<<endl;
        lrs_->randomizeGids();

        lrs_->fillWithZeroCenters(ct.numst);
        
        if( onpe0 && ct.verbose>0 )os_<<"readInput: Read "<<count
                      <<" orbital centers from input file"<<endl;
    } // mmpi.instancePE0()
   
    if( ct.verbose>0 )
        printWithTimeStamp("setup LRs...",os_);
    lrs_->setup();
    lrs_->printInfo(os_);

#ifdef DEBUG
    lrs_->print(os_);
#endif
    
    return ct.numst;
}

int MGmol::readCoordinates(ifstream* tfile,
                           const bool cell_relative)
{     
    Control& ct = *(Control::instance());
    if( ct.verbose>0 )
        printWithTimeStamp("Read atomic coordinates...",os_);
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();

    const double lattice[3]={mygrid.ll(0),mygrid.ll(1),mygrid.ll(2)};

    // setup ions
    const std::vector<Species>& sp(ct.getSpecies());
    ions_=new Ions(lattice,sp);

    if( ct.restart_info>0 && ct.override_restart==0 ) //read restart ionic positions
    {
        if( ct.restart_info>0 )
        {
            if(onpe0 && ct.verbose>0)
            {
                 os_<<"Initialize ionic positions from restart file "
                    <<ct.restart_file<<endl;
            }
            ions_->initFromRestartFile(*h5f_file_);
        }   
    }
    else
    {
        // Coordinates and species type for each ion. 
        int info=ions_->readAtoms(tfile,cell_relative);
        
        return info;
    }
    
    const int num_ions=ions_->getNumIons();
    if(onpe0)os_<<num_ions<<" ions in simulation"<<endl;

    return 0;
}

int MGmol::readCoordinates(const string filename,
                           const bool cell_relative)
{ 
    Control& ct = *(Control::instance());
    if( ct.verbose>0 )
        printWithTimeStamp("Read atomic coordinates...",os_);
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();

    const double lattice[3]={mygrid.ll(0),mygrid.ll(1),mygrid.ll(2)};

    // setup ions
    const std::vector<Species>& sp(ct.getSpecies());
    ions_=new Ions(lattice,sp);

    if( ct.restart_info>0 && ct.override_restart==0 ) //read restart ionic positions
    {
        if( ct.restart_info>0 )
        {
            if(onpe0 && ct.verbose>0)
            {
                 os_<<"Initialize ionic positions from restart file "
                    <<ct.restart_file<<endl;
            }
            ions_->initFromRestartFile(*h5f_file_);
        }   
    }
    else
    {
        // Coordinates and species type for each ion. 
        int info=ions_->readAtoms(filename,cell_relative);
        
        return info;
    }
    
    const int num_ions=ions_->getNumIons();
    if(onpe0)os_<<num_ions<<" ions in simulation"<<endl;

    return 0;
}
    
// Reads and parses the input control file 
int MGmol::readInput(const string input_file)
{
    printWithTimeStamp("MGmol::readInput()...",cout);

    ifstream* tfile=0;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Control& ct = *(Control::instance());

    if( mmpi.instancePE0() ){
        if( onpe0 )os_<<"Read input file "<<input_file<<endl;
        tfile = new ifstream(input_file.data(), ios::in);
        if ( !tfile->is_open() )
        {
            cerr << " Unable to open file " << input_file.data() << endl;
            global_exit(0);
        }else{
            if( ct.verbose>0 )
                os_<<"Open "<<input_file.data()<<endl;
        }
    }
    
    bool cell_relative=false;
    int status=readParameters(tfile, cell_relative);
    if( status==-1 )return -1;

    Mesh* mymesh = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    if( ct.restart_info>0 )
        h5f_file_=new HDFrestart(ct.restart_file, myPEenv, ct.restart_file_type);

    status=readCoordinates(tfile, cell_relative);
    if( status==-1 )return -1;
    
    const short myspin=mmpi.myspin();
    const int nval=ions_->getNValenceElectrons();
    ct.setNumst(myspin,nval);
    ct.setTolEnergy();
    
    // create localization regions
    const pb::Grid& mygrid  = mymesh->grid();
    Vector3D vcell(mygrid.ll(0),mygrid.ll(1),mygrid.ll(2));
    lrs_=new LocalizationRegions(vcell,ct.tol_orb_centers_move);

    readLRsFromInput(tfile);

    constraints_->readConstraints(tfile);

    if( mmpi.instancePE0() )
    {
        os_<<"Close "<<input_file.data()<<endl;
        tfile->close();
        delete tfile;
    }
    
    /* define spread radius - this is done last 
     * to ensure that required data is available 
    */
    ct.setSpreadRadius();

    return 0;
}


    

// Reads and parses the input control file 
int MGmol::readInput(const string input_file1,const string input_file2)
{
    ifstream* tfile1=0;
    ifstream* tfile2=0;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Mesh* mymesh = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    Control& ct = *(Control::instance());

    printWithTimeStamp("MGmol::readInput()...",os_);
    if(mmpi.instancePE0()){
        os_<<"Read input file "<<input_file1<<endl;
        tfile1 = new ifstream(input_file1.data(), ios::in);
        if ( !tfile1->is_open() )
        {
            cerr << " Unable to open file " << input_file1.data() << endl;
            global_exit(0);
        }else{
            os_<<"Open "<<input_file1.data()<<endl;
        }
    }
    if(mmpi.instancePE0()){
        os_<<"Read input file "<<input_file2<<endl;
        tfile2 = new ifstream(input_file2.data(), ios::in);
        if ( !tfile2->is_open() )
        {
            cerr << " Unable to open file " << input_file2.data() << endl;
            global_exit(0);
        }else{
            os_<<"Open "<<input_file2.data()<<endl;
        }
    }

    bool cell_relative=false;
    readParameters(tfile1, cell_relative);

    if( ct.restart_info>0 )
        h5f_file_=new HDFrestart(ct.restart_file, myPEenv, ct.restart_file_type);

    readCoordinates(tfile2, cell_relative);

    const short myspin=mmpi.myspin();
    const int nval=ions_->getNValenceElectrons();
    ct.setNumst(myspin,nval);
    
    // create localization regions
    const pb::Grid& mygrid  = mymesh->grid();
    Vector3D vcell(mygrid.ll(0),mygrid.ll(1),mygrid.ll(2));
    lrs_=new LocalizationRegions(vcell,ct.tol_orb_centers_move);

    readLRsFromInput(tfile2);

    constraints_->readConstraints(tfile1);

    if( mmpi.instancePE0() )
    {
        os_<<"Close "<<input_file1.data()<<endl;
        os_<<"Close "<<input_file2.data()<<endl;

        tfile1->close();
        tfile2->close();
        delete tfile1,tfile2;
    }

    /* define spread radius if needed - this is done last 
     * to ensure that required data is available 
    */
    ct.setSpreadRadius();
    printWithTimeStamp("MGmol::readInput()... done...",os_);

    return 0;
}


