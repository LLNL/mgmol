// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
using namespace std;

#include "global.h"

#include "MGmol.h"
#include "Forces.h"
#include "GridFunc.h"
#include "FDoper.h"
#include "LocGridOrbitals.h"
#include "MLWFTransform.h"
#include "Ions.h"
#include "MPIdata.h"
#include "SparseDistMatrix.h"
#include "SubMatrices.h"
#include "Hamiltonian.h"
#include "Control.h"
#include "Rho.h"
#include "XConGrid.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "AndersonMix.h"
#include "ProjectedMatricesMehrstellen.h"
#include "ProjectedMatricesSparse.h"
#include "HDFrestart.h"
#include "ConstraintSet.h"
#include "Mesh.h"
#include "Potentials.h"
#include "LocalizationRegions.h"
#include "PBEonGrid.h"
#include "PBEonGridSpin.h"
#include "LDAonGrid.h"
#include "KBPsiMatrix.h"
#include "KBPsiMatrixSparse.h"
#include "MasksSet.h"
#include "DistMatrix.h"
#include "LBFGS.h"
#include "FIRE.h"
#include "ABPG.h"
#include "GrassmanLineMinimization.h"
#include "SubspaceProjector.h"
#include "SpreadsAndCenters.h"
#include "SpreadPenalty.h"
#include "EnergySpreadPenalty.h"
#include "SpreadPenaltyVolume.h"
#include "AOMMprojector.h"
#include "DFTsolver.h"
#include "DMStrategyFactory.h"
#include "Preconditioning.h"
#include "Masks4Orbitals.h"
#include "OrbitalsPreconditioning.h"
#include "PackedCommunicationBuffer.h"
#include "PoissonInterface.h"
#include "MDfiles.h"

#define DELTA  1e-8

namespace mgmol{
    std::ostream* out = NULL;
}

dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE> *MGmol::remote_tasks_DistMatrix_ptr_=0;

string description;

extern Timer computeHij_tm;
extern Timer get_kbpsi_tm;
extern Timer get_Hpsi_and_Hij_tm;
extern Timer vnlpsi_tm;
extern Timer sygv_tm;
extern Timer split_allreduce_sums_double_tm;
extern Timer sgemm_tm;
extern Timer dgemm_tm;
extern Timer mpgemm_tm;
extern Timer tttgemm_tm;
extern Timer dsyrk_tm;
extern Timer ssyrk_tm;
extern Timer mpsyrk_tm;
extern Timer tttsyrk_tm;
extern Timer mpdot_tm;
extern Timer ttdot_tm;
extern Timer get_NOLMO_tm;
extern Timer get_MLWF_tm;
extern Timer md_iterations_tm;
extern Timer md_tau_tm;
extern Timer md_moveVnuc_tm; 
extern Timer md_updateMasks_tm;        
extern Timer md_extrapolateOrbitals_tm;         
extern Timer md_updateRhoAndPot_tm;
extern Timer quench_tm;
extern Timer ions_setupInteractingIons_tm;
extern Timer ions_setup_tm;
extern Timer updateCenters_tm;
//extern Timer nonOrthoRhoKernel_tm;
//extern Timer nonOrthoRhoKernelDiagonalBlock_tm;

Timer MGmol::total_tm_("MGmol::total");
Timer MGmol::setup_tm_("MGmol::setup");
Timer MGmol::closing_tm_("MGmol::closing");
Timer MGmol::init_tm_("MGmol::init");
Timer MGmol::dump_tm_("MGmol::dump");
Timer MGmol::evnl_tm_("MGmol::evnl");
Timer MGmol::get_res_tm_("MGmol::comp_res_from_Hphi");
Timer MGmol::computeHij_tm_("MGmol::computeHij");
Timer MGmol::get_Hpsi_and_Hij_tm_("MGmol::get_Hpsi_and_Hij");
Timer MGmol::comp_res_tm_("MGmol::comp_res");
Timer MGmol::init_nuc_tm_("MGmol::init_nuc");

#include "Signal.h"
set<int> Signal::recv_;

MGmol::MGmol(MPI_Comm comm, std::ostream& os):
    os_(os)
{
    comm_=comm;
    
    constraints_=new ConstraintSet();
    
    mgmol::out = &os;

    geom_optimizer_=0;
    lrs_=0;
    local_cluster_=0;
    proj_matrices_=0;
    dm_strategy_=0;
    
    h5f_file_=0;
    
    aomm_=0;

    spreadf_=0;
    
    spread_penalty_=0;
    
    data_distributor_ = new BasicDataDistributors();
 
    orbitals_precond_=0;
    
    forces_=0;
    
    energy_=0;
}
    
MGmol::~MGmol()
{
    Control& ct = *(Control::instance());
    delete[] rhoc_;
    delete electrostat_;
    delete rho_;
    delete constraints_;
    if( !ct.short_sighted )
    {
      delete remote_tasks_DistMatrix_;
      remote_tasks_DistMatrix_=0;
    }
    delete xcongrid_;
    delete energy_;
    if( hamiltonian_!=0 )delete hamiltonian_;
    if( geom_optimizer_!=0 )delete geom_optimizer_;
    
    delete currentMasks_;
    delete corrMasks_;

    if( aomm_!=0 )delete aomm_;

    delete current_orbitals_;
    delete ions_;
    delete g_kbpsi_;
    
    delete proj_matrices_;
    delete lrs_;
    if(local_cluster_!=0) delete local_cluster_;
    
    if( h5f_file_!=0 )delete h5f_file_;
    
    if( spreadf_!=0 )delete spreadf_;
    
    if( spread_penalty_!=0 )delete spread_penalty_;
    
    delete data_distributor_;
    
    delete forces_;
    if(dm_strategy_!=0) delete dm_strategy_;
}


int MGmol::initial()
{  
    Control& ct = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Mesh* mymesh = Mesh::instance();
    
    assert( ct.numst>=0 );
    
    if( ct.verbose>0 )
        printWithTimeStamp("MGmol::initial()...",os_);

    init_tm_.start();
    
    const pb::Grid& mygrid  = mymesh->grid();

    hamiltonian_->setup(mygrid,ct.lap_type);
    
    pb::Lap<ORBDTYPE>* lapop = ct.Mehrstellen() ? hamiltonian_->lapOper() : 0;
    g_kbpsi_ = new KBPsiMatrixSparse(lapop);
    
    check_anisotropy();
    
    if( ct.verbose>0 )
        printWithTimeStamp("MGmol::initial(), create ProjectedMatrices...",os_);

    // If not an initial run read data from files  
    if( ct.restart_info>2 )
    {
        string name = "ExtrapolatedFunction";
        if( ct.verbose>0 )
        printWithTimeStamp("read LRs from ExtrapolatedFunction database...",os_);
        int n=read_restart_lrs(*h5f_file_, name);
        if( n==0 )
        {
            if( ct.verbose>0 )
            printWithTimeStamp("read LRs from Function database...",os_);
            name="Function";
            n=read_restart_lrs(*h5f_file_, name);
        }
        if( n<0 )return n;
        
        if( n>0 )lrs_->setup();
    }

    double dlrsmin=lrs_->computeMinDistBetweenLocalPairs();
    if(dlrsmin<1.e-3)
    {
        cerr<<"ERROR: Min. distance between LR centers is "
            <<dlrsmin<<"!!!"<<endl;
        return 1;
    }

    // initialize and setup load balancing object
    // set number of iterations to 10.
    if(ct.load_balancing_alpha > 0.0)
    {
       local_cluster_ = new ClusterOrbitals(lrs_);
       local_cluster_->setup();
       local_cluster_->computeClusters(ct.load_balancing_max_iterations);
    }
    // initialize data distribution objects
    const pb::PEenv& myPEenv=mymesh->peenv();
    double domain[3]={mygrid.ll(0),mygrid.ll(1),mygrid.ll(2)};
    data_distributor_->initialize(lrs_, myPEenv, domain);   
    
    bool with_spin=(mmpi.nspin()>1);
    if( ct.Mehrstellen() )
        proj_matrices_ = new ProjectedMatricesMehrstellen(ct.numst, with_spin );
    else if(ct.short_sighted)
        proj_matrices_ = new ProjectedMatricesSparse(ct.numst, lrs_, local_cluster_); 
    else
        proj_matrices_ = new ProjectedMatrices(ct.numst, with_spin );

    forces_ = new Forces(hamiltonian_,rho_,proj_matrices_);
    
    if( ct.verbose>0 )
        printWithTimeStamp("MGmol::initial(), create MasksSet...",os_);

    currentMasks_=new MasksSet(false,ct.getMGlevels());
    currentMasks_->setup(*lrs_);

    corrMasks_=new MasksSet(true,0);
    corrMasks_->setup(*lrs_);

    BlockVector<ORBDTYPE>::setOverAllocateFactor(
        ct.orbitalsOverallocateFactor() );

    if( ct.verbose>0 )
        printWithTimeStamp("MGmol::initial(), create LocGridOrbitals...",os_);

    current_orbitals_=new LocGridOrbitals(mygrid, mymesh->subdivx(),
                                         ct.numst,
                                         ct.bc,
                                         proj_matrices_,
                                         lrs_,
                                         currentMasks_,
                                         corrMasks_,
                                         local_cluster_);    
    current_orbitals_->setup();    
    
    if( !ct.short_sighted )
    {
        printWithTimeStamp("MGmol::initial(), create MatricesBlacsContext and misc...",os_);

        MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
        const dist_matrix::BlacsContext& bc=*mbc.bcxt();
        dist_matrix::DistMatrix<DISTMATDTYPE> tmp("tmp", bc, ct.numst, ct.numst);
        remote_tasks_DistMatrix_
            =new dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>(tmp);
        ProjectedMatrices::registerRemoteTasksDistMatrix(remote_tasks_DistMatrix_);
        LocalMatrices<MATDTYPE>::registerRemoteTasksDistMatrix(remote_tasks_DistMatrix_);
//        g_kbpsi_->registerRemoteTasksDistMatrix(remote_tasks_DistMatrix_);
        remote_tasks_DistMatrix_ptr_ = remote_tasks_DistMatrix_;
    }
    if( ct.wf_extrapolation==0 )
    {
        BlockVector<ORBDTYPE>::incMaxAllocInstances(2);
    }    
    if( ct.it_algo_type==3 )
    {
        BlockVector<ORBDTYPE>::incMaxAllocInstances(2);
    }    
    if( ct.wf_m>1 || ct.wf_extrapolation>1 )
        BlockVector<ORBDTYPE>::incMaxAllocInstances(1);
    for(short i=1;i<ct.wf_m; i++)
        BlockVector<ORBDTYPE>::incMaxAllocInstances(2);
    if( ct.use_kernel_functions )
        BlockVector<ORBDTYPE>::incMaxAllocInstances(1);
    
    Potentials& pot =hamiltonian_->potential();
    pb::Lap<ORBDTYPE>* lapOper=hamiltonian_->lapOper();

    // If not an initial run read data from files  
    if( ct.restart_info>0 )
    {
        md_time_     =h5f_file_->getMDTimeFromFile();
        md_iteration_=h5f_file_->getMDstepFromFile();
        
        int ierr=read_restart_data(*h5f_file_, *rho_,*current_orbitals_);
        if( ierr<0 )return -1;

        if( ct.restart_info>2 )
        {
            current_orbitals_->applyMask();
        }
    }
    else
    {
        md_time_=0.;
        md_iteration_=1;
    }
    
    mmpi.barrier();
    if( ct.verbose>0 )current_orbitals_->printChromaticNumber(os_);
    
    initBackground();

    // Random initialization of the wavefunctions
    if( ct.restart_info<=2 )
    {
        if( ct.verbose>0 ) printWithTimeStamp("MGmol::initial(), init wf and masks...",os_);

        // Make temp mask for initial random wave functions  
        if( ct.init_loc==1 )
        {
            float  cut_init=ct.initRadius();
            assert(cut_init>0.);
            currentMasks_->initialize(*lrs_, 0, cut_init);
        }
    
        // Initialize states
        current_orbitals_->initWF();
        
        // initialize masks again
        if( ct.init_loc==1 )
        {
            currentMasks_->initialize(*lrs_, 0);
        }
    }

    // Initialize the radial potential stuff  
    if( ct.verbose>0 )
        printWithTimeStamp("initKBR()...",os_);
    initKBR();

    ions_->setup();

    double d=ions_->computeMinLocalSpacing();
    if( d<1.e-3 )
    {
       cerr<<"ERROR: min. distance between centers is smaller than 1.e-3!!!\n";
    }
    
    // Initialize the nuclear local potential and the compensating charges  
    if( ct.verbose>0 )
        printWithTimeStamp("initNuc()...",os_);
    initNuc(*ions_);

    // initialize Rho
    if( ct.verbose>0 )
        printWithTimeStamp("Initialize Rho...",os_);
    if( ct.restart_info <=1 )
        proj_matrices_->setDMuniform(ct.getNelSpin(),current_orbitals_->getIterativeIndex());
    
    rho_->setup(ct.orbital_type, current_orbitals_->getGlobalIndexes());
    
    if( ct.restart_info <=1 )
    {
        if( ct.init_type==0 ) // random functions
        {
            rho_->init(rhoc_);
        }
        else
        {
            rho_->update(*current_orbitals_);
        }
    }
    
    if( ct.verbose>0 )
        printWithTimeStamp("Initialize XC functional...",os_);
    if( ct.xctype==0 )
    {
        xcongrid_=new LDAonGrid(*rho_,pot);
    }
    else if( ct.xctype==2 )
    {
        if( mmpi.nspin()>1 )
            xcongrid_=new PBEonGridSpin(*rho_,pot);
        else
            xcongrid_=new PBEonGrid(*rho_,pot);
    }
    else
    {
        (*MPIdata::serr)<<"Invalid XC option"<<endl;
        return -1;
    }

    // initialize nl potentials with restart values if possible
//    if( ct.restart_info>1 )
    {
        electrostat_->setupInitialVh(pot.vh_rho());
        electrostat_->computeVhRho(*rho_); // for energy
        
        // xc computed from rho
        xcongrid_->update();
    }

    current_orbitals_->setDataWithGhosts();
    current_orbitals_->trade_boundaries();

//    if(ct.restart_info <= 1)pot.initWithVnuc();
    // initialize matrices S and invB      

    if( ct.numst>0 )
    {
        if( ct.verbose>0 )
            printWithTimeStamp("Compute initial Gram matrix...",os_);
        current_orbitals_->computeGramAndInvS(ct.verbose);
        if( ct.verbose>0 )
            printWithTimeStamp("Compute initial B matrix...",os_);
        current_orbitals_->computeBAndInvB(*lapOper);        
        if( ct.verbose>0 )
            printWithTimeStamp("Compute initial condition number...",os_);
        current_orbitals_->checkCond(100000., (ct.atoms_dyn!=0));
    }
        
    if( ct.verbose>0 )
        printWithTimeStamp("Setup kbpsi...",os_);
    g_kbpsi_->setup(*ions_, *current_orbitals_);
    
    if( ct.restart_info == 0 )
    {
        if( ct.verbose>0 )
            printWithTimeStamp("update_pot...",os_);
        update_pot(*ions_);
    }

    if( ct.isLocMode() || ct.isSpreadFunctionalActive() )
    {
        Vector3D origin(mygrid.origin(0),
                        mygrid.origin(1),
                        mygrid.origin(2));
        Vector3D ll(mygrid.ll(0),mygrid.ll(1),mygrid.ll(2));
        spreadf_ = new SpreadsAndCenters(origin,ll);
    }
    
    bool energy_with_spread_penalty=false;
    if( ct.isSpreadFunctionalActive() )
    {
        if( ct.isSpreadFunctionalVolume() )
        {
            spread_penalty_=new SpreadPenaltyVolume(
                spreadf_,
                ct.spreadPenaltyTarget(),
                ct.spreadPenaltyAlphaFactor(),
                ct.spreadPenaltyDampingFactor());
        
        }
        else
        if( ct.isSpreadFunctionalEnergy() )
        {
            energy_with_spread_penalty=true;
            spread_penalty_=new EnergySpreadPenalty(
                spreadf_,
                ct.spreadPenaltyTarget(),
                ct.spreadPenaltyAlphaFactor());
        
        }
        else
            spread_penalty_=new SpreadPenalty(
                spreadf_,
                ct.spreadPenaltyTarget(),
                ct.spreadPenaltyAlphaFactor(),
                ct.spreadPenaltyDampingFactor());
    }
    
    SpreadPenaltyInterface* spread_penalty = energy_with_spread_penalty ? spread_penalty_ : 0;
    energy_=new Energy(mygrid,*ions_,pot,*electrostat_,*rho_,*xcongrid_,spread_penalty);

    if( ct.verbose>0 )
        printWithTimeStamp("Setup matrices...",os_);
    
    updateHmatrix(*current_orbitals_, *ions_);
    
    // HMVP algorithm requires that H is initialized
    dm_strategy_ = DMStrategyFactory::create(comm_, os_, 
        *ions_,
        rho_,
        energy_,
        electrostat_,
        this,
        proj_matrices_,current_orbitals_);
    


    // theta = invB * Hij 
    proj_matrices_->updateThetaAndHB();

    dm_strategy_->initialize();

    if( ct.verbose>1 )
    {
        proj_matrices_->printMatrices(os_);
        proj_matrices_->printDM(os_);
    }

    init_tm_.stop();

    if( ct.verbose>0 )
        printWithTimeStamp("Initialization done...",os_);

    return 0;
} // initial()

void MGmol::run()
{
    total_tm_.start();

    setup();
      
    Control& ct = *(Control::instance());

    double eks=0.;
    
    if( ct.verbose>0 )
        printWithTimeStamp("Run...",os_);
    // Dispatch to the method chosen 
    switch(ct.atoms_dyn)
    {
        case 0:                  // Quench the electrons
            quench(current_orbitals_,
                   *ions_, ct.max_electronic_steps, 20, eks);

            // Forces for the last states 
            force(*current_orbitals_, *ions_);
            
            constraints_->addConstraints(*ions_);

            constraints_->setup(*ions_);
            
            constraints_->printConstraintsForces(os_);

            constraints_->projectOutForces(20);
            
            if( (ions_->getNumIons()<=1024 || ct.verbose>2) && ct.verbose>0 )
                ions_->printForcesGlobal(os_);

            finalEnergy();
            
            printMM();
            
            break; 
    
        case 2:                // MD
            md(&current_orbitals_, *ions_);
            //finalEnergy();
            break;
        
        case 6:                // LBFGS
            lbfgsrlx(&current_orbitals_,*ions_);
            //finalEnergy();
            break;

        case 7:                // FIRE
            runfire(&current_orbitals_,*ions_);
            //finalEnergy();
            break;

        default:
            (*MPIdata::serr)<<"run: Undefined MD method"<<endl;

    }  

    cleanup();
} //run()
  

void MGmol::finalEnergy()
{
    // Get the total energy
    const double ts=0.5*proj_matrices_->computeEntropy(); // in [Ha]
    total_energy_ = energy_->evaluateTotal(ts, proj_matrices_, *current_orbitals_, 2, os_);
}
    
void MGmol::printMM()
{
    Control& ct = *(Control::instance());
    if( ct.tmatrices==1 )
    {
        ofstream tfile("s.mm", ios::out);
        proj_matrices_->printGramMM(tfile);
        ofstream tfileh("h.mm", ios::out);
        ProjectedMatrices* projmatrices=dynamic_cast<ProjectedMatrices*>(proj_matrices_);
        assert(projmatrices!=0);
        projmatrices->printHamiltonianMM(tfileh);
    }
}

/* Writes out header information */
void MGmol::write_header()
{
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();

    time_t tt;
  
    time( &tt );
    char* timeptr = ctime( &tt );

    Control& ct = *(Control::instance());

    if( onpe0 )
    {
    
    os_<<"//////////////////////////////////////////////////////////"<<endl;
    os_<<endl;
#ifdef GITHASH
#define xstr(x) #x
#define LOGGIT(x) os_<<" MGmol: git_hash "<<xstr(x)<<endl;
    LOGGIT(GITHASH);
    os_<<endl;
#endif
    os_<<" Compiled: "<<__DATE__<<", "<<__TIME__<<endl;
    os_<<" Real-space finite difference ab initio calculations"<<endl;
    os_<<endl;
    os_<<" authors: J.-L. Fattebert, Lawrence Livermore National Laboratory"<<endl;
    os_<<"          D. Osei-Kuffuor, Lawrence Livermore National Laboratory"<<endl;
    os_<<endl;
    os_<<"//////////////////////////////////////////////////////////"<<endl;

    os_<<endl<<endl<<description<<endl;
    os_<<" Run started at "<<timeptr<<endl;
    
    os_<<" Orbitals precision (in bytes): "<<sizeof(ORBDTYPE)<<endl<<endl;

    Potentials& pot =hamiltonian_->potential();
    pot.writeNames(os_);

    mymesh->print(os_);
    
    pb::Lap<ORBDTYPE>* lapop=hamiltonian_->lapOper();
    os_<<" Laplacian Discretization: "<<lapop->name()<<endl;

#ifdef SMP_NODE
    os_ << " " << omp_get_max_threads() << " thread"
         << ( omp_get_max_threads() > 1 ? "s " : " " );
    os_ << "active" << endl << endl;
#endif
    
    os_<<" ScaLapack block size: "<<dist_matrix::DistMatrix<DISTMATDTYPE>::getBlockSize()<<endl;
    
    if( !ct.short_sighted )
    {
        MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
        mbc.print(os_);
    }
    } // onpe0
    
    if( current_orbitals_!=NULL && ct.verbose>0 )
    {
        current_orbitals_->printNumst(os_);
        current_orbitals_->printChromaticNumber(os_);
    }
    
    int nions=ions_->getNumIons();
    if( onpe0 )
    {
    os_<<" Number of ions     = "<<nions<<endl;
    os_<<" Total charge in cell = "<<charge_in_cell_<<endl;
 
    ct.print(os_);
    
    os_<<endl<<endl<<" Atomic species information:"<<endl<<endl;
    os_<<fixed<<setprecision(5);
    
    const std::vector<Species>& sp(ct.getSpecies());
    for(int idx = 0;idx < (int)sp.size();idx++)
    {
        const Species& sp_ion(sp[idx]);
    
        os_<<" Species #"<<idx + 1<<": ";
        sp_ion.print(os_);
 
        os_<<" dim_l        = "<<sp_ion.dim_l()
            <<"    ->Diameter    = "<<sp_ion.dim_l()*mygrid.hmin()<<"[bohr]"<<endl;
        os_<<" dim_nl      = "<<sp_ion.dim_nl()
            <<"    ->Diameter    = "<<sp_ion.dim_nl()*mygrid.hmin()<<"[bohr]"<<endl;
 
    }    
    if(ct.short_sighted)
    {  
       os_<<endl;
       os_<<" Short_Sighted Solver Parameters: "<<endl<<endl;
       os_<<" Accelerator = Flexible GMRES("<<ct.fgmres_kim<<")"<<endl;
       if(ct.ilu_type==0)  os_<<" Preconditioner = ilu"<<ct.ilu_lof<<endl;
       else if(ct.ilu_type==1)  os_<<" Preconditioner = modified ilu"<<endl;       
       else if(ct.ilu_type==2)  os_<<" Preconditioner = standard ilu"<<endl;              
       else os_<<" Preconditioner = ilu0 preconditioner"<<endl;
       os_<<scientific<<" fgmres convergence tolerance = "<<ct.fgmres_tol<<endl;
       os_<<" fgmres max. iterations = "<<ct.fgmres_maxits<<endl;
       os_<<scientific<<" drop tolerance for (standard/ modified) ilu = "<<ct.ilu_droptol<<endl;
       os_<<" max. fill-in for (standard/ modified) ilu = "<<ct.ilu_maxfil<<endl;     
    } 
  }// onpe0

  // Write out the ionic postions and displacements  
  ions_->printPositions(os_);   
    
  if( current_orbitals_!=NULL && ct.verbose>3 )
      lrs_->printAllRegions(os_);

} 

void MGmol::global_exit(int i)
{
#ifdef USE_MPI
    MPI_Abort(comm_, i);
#else
    abort();
#endif
}


void MGmol::check_anisotropy()
{
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();

    if(mygrid.anisotropy() > 1.2)
    {
        write_header();
        if( onpe0 )
            (*MPIdata::serr)<<" hmax="<<mygrid.hmax()
                <<", hmin="<<mygrid.hmin()<<endl;
        (*MPIdata::serr)<<"init: Anisotropy too large: "<<mygrid.anisotropy()<<endl;
        global_exit(2);
    }

}

void MGmol::initBackground()
{
    assert( ions_!=0 );
    
    Control& ct = *(Control::instance());

    // Count up the total ionic charge  
    ionic_charge_ = ions_->computeIonicCharge();

    // calculation the compensating background charge
    //   for charged supercell calculations  
    background_charge_ = 0.;
    charge_in_cell_ = ionic_charge_ - ct.getNel();
    if( ct.bcPoisson[0]!=2 && ct.bcPoisson[1]!=2 && ct.bcPoisson[2]!=2 )
    {
        background_charge_ = (-1.)*charge_in_cell_;
    }
    if( onpe0 && ct.verbose>0 )
    {
        os_<<"N electrons=      "<<ct.getNel()<<endl;
        os_<<"ionic charge=     "<<ionic_charge_<<endl;
        os_<<"background charge="<<background_charge_<<endl;
    }

    if(fabs(background_charge_)<1.e-10)background_charge_=0.;
}

void MGmol::printEigAndOcc()
{
    Control& ct = *(Control::instance());
    if( ! (ct.fullyOccupied() 
        && ct.orbital_type>0 
        && ct.occupationWidthIsZero())
        && onpe0 )
    {
        proj_matrices_->printEigenvalues(os_);
        proj_matrices_->printOccupations(os_);
    }
}

double MGmol::get_charge(RHODTYPE *rho)
{
    Control& ct = *(Control::instance());
    Mesh* mymesh = Mesh::instance();
    const pb::PEenv& myPEenv=mymesh->peenv();
    const pb::Grid& mygrid  = mymesh->grid();

    double  charge = 0.;
    const int numpt=(int)mygrid.size();
    
    for(int idx = 0;idx < numpt;idx++)
    {
        assert( rho[idx]<100. );
        assert( rho[idx]>=-1. ); // coulb be < 0. with neutralizing background
        charge += rho[idx];
    } 
    charge = myPEenv.double_sum_all(charge);

    assert(mygrid.vel()>0.000001);
    assert(mygrid.vel()<1000.);

    charge *= mygrid.vel();    
    if( onpe0 && ct.verbose>0 )
        os_<<setprecision(8)<<fixed<<"Charge: "<<charge<<endl;
    
    return charge;
}

//  Initialization of the compensating potential
void MGmol::initVcomp(Ions& ions)
{
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();

    double  point[3]={0.,0.,0.};
    Control& ct = *(Control::instance());

    Potentials& pot =hamiltonian_->potential();
    pot.set_vcomp(0.);

    const int dim0=mygrid.dim(0);
    const int dim1=mygrid.dim(1);
    const int dim2=mygrid.dim(2);

    const double start0=mygrid.start(0);
    const double start1=mygrid.start(1);
    const double start2=mygrid.start(2);

    const double h0=mygrid.hgrid(0);
    const double h1=mygrid.hgrid(1);
    const double h2=mygrid.hgrid(2);
    
    const double lattice[3]={mygrid.ll(0),mygrid.ll(1),mygrid.ll(2)};
    
    const int incx=dim2*dim1;

    /* Loop over ions */
    vector<Ion*>::const_iterator ion=ions.overlappingVL_ions().begin();
    while(ion!=ions.overlappingVL_ions().end())
    {
        const double lrad = (*ion)->getRadiusLocalPot();
        assert( lrad>0.1 );
    
        const double Zv = (*ion)->getZion();
        const double rc = (*ion)->getRC();
        assert( rc>0.1 );
        const double invrc = 1./rc;

        for(int ix = 0;ix < dim0;ix++)
        {
            point[0] = start0+ix*h0;
            int istart = incx*ix;

            for(int iy = 0;iy < dim1;iy++)
            {
                point[1] = start1+iy*h1;
                int jstart = istart + dim2*iy;

                for(int iz = 0;iz < dim2;iz++)
                {
                    point[2] = start2+iz*h2;

                    const double r = (*ion)->minimage(point,lattice, ct.bcPoisson);

                    if(r< lrad)
                    {
                        if(r <= DELTA )
                        {
                            pot.add_vcomp(jstart+iz, Zv * M_2_SQRTPI * invrc);
                        }
                        else
                        {
                            pot.add_vcomp(jstart+iz, Zv * erf(r*invrc) / r);
                        }
                    }
 
                } 

            } 

        } 

        ion++;
    }
}

#if 0
double get_trilinval(const double xc, const double yc, const double zc,
                     const double h0, const double h1, const double h2,
                     const Vector3D& ref, const Vector3D& lattice,
                     RadialInter& lpot)
{
    double val=0.;
    const double shift=0.5;
    const double hh0=shift*h0;
    const double hh1=shift*h1;
    const double hh2=shift*h2;
    for(int i=-1;i<2;i++){
        double xcc=xc+i*hh0;
        double w0=1.-fabs(0.5*i);
        for(int j=-1;j<2;j++){
            double ycc=yc+j*hh1;
            double w1=1.-fabs(0.5*j);
            for(int k=-1;k<2;k++){
                double zcc=zc+k*hh2;
                double w2=1.-fabs(0.5*k);
                
                Vector3D vc(xcc,ycc,zcc);
                double w=w0*w1*w2*0.125;
                double r = vc.minimage(ref,lattice);            
                val += w*lpot.cubint(r);
            }
        }
    }
    return val;
}
#endif


void MGmol::initNuc(Ions& ions)
{
    init_nuc_tm_.start();
    
    Control& ct = *(Control::instance());
    if( onpe0 && ct.verbose>1 )
        os_<<"Init vnuc and rhoc..."<<endl;
    
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();

    Potentials& pot =hamiltonian_->potential();
    POTDTYPE* vnuc=pot.vnuc();
    
    const int numpt=mygrid.size();

    assert(numpt>0);
    assert(mygrid.vel()>0.);

    memset(rhoc_, 0, numpt*sizeof(RHODTYPE));
    memset(vnuc, 0, numpt*sizeof(POTDTYPE));
    
    const double h0=mygrid.hgrid(0);
    const double h1=mygrid.hgrid(1);
    const double h2=mygrid.hgrid(2);
    
    Vector3D ll(mygrid.ll(0),mygrid.ll(1),mygrid.ll(2));
    
    const int dim0=mygrid.dim(0);
    const int dim1=mygrid.dim(1);
    const int dim2=mygrid.dim(2);
    const int gdim0=mygrid.gdim(0);
    const int gdim1=mygrid.gdim(1);
    const int gdim2=mygrid.gdim(2);
    
    const int inc0=dim1*dim2;
    
    const double pi3half=M_PI*sqrt(M_PI);

    const int ilow = mygrid.istart(0);
    const int jlow = mygrid.istart(1);
    const int klow = mygrid.istart(2);
    const int ihi = ilow + dim0 - 1;
    const int jhi = jlow + dim1 - 1;
    const int khi = klow + dim2 - 1;

    // Loop over ions
    vector<Ion*>::const_iterator ion=ions.overlappingVL_ions().begin();
    while(ion!=ions.overlappingVL_ions().end())
    {
        assert( (*ion)->map_l() );
        
        Vector3D position((*ion)->position(0),(*ion)->position(1),(*ion)->position(2));

        // Generate indices
        vector<vector<int> >  Ai;
        Ai.resize(3);
        (*ion)->get_Ai(Ai[0],gdim0,0);
        (*ion)->get_Ai(Ai[1],gdim1,1);
        (*ion)->get_Ai(Ai[2],gdim2,2);
        const int dimlx=Ai[0].size();
        const int dimly=Ai[1].size();
        const int dimlz=Ai[2].size();

        const double  Zv = (*ion)->getZion();
        const double  rc = (*ion)->getRC();
        const double  lr = (*ion)->getRadiusLocalPot();
        assert(Zv>=0.);
        assert(rc>1.e-8);

        const double  rcnorm = rc*rc*rc*pi3half;
        assert(rcnorm>1.e-8);
        const double  alpha=Zv/rcnorm;
        const double  inv_rc2=1./(rc * rc);

        const RadialInter& lpot( (*ion)->getLocalPot() );
        int icount = 0;
        double charge_ion=0.;
        double zc = (*ion)->lstart(2);
    
        for(int iz = 0;iz <  dimlz;iz++)
        {
            double yc = (*ion)->lstart(1);
            const int ai2iz=Ai[2][iz];

            if( (ai2iz >= klow) && (ai2iz <= khi) )
            for(int iy = 0;iy <  dimly;iy++)
            {
                double xc = (*ion)->lstart(0);
                const int ai1iy=Ai[1][iy];

                if( (ai1iy >= jlow) && (ai1iy <= jhi) )
                for(int ix = 0;ix <  dimlx;ix++)
                {
                    const int ai0ix=Ai[0][ix];
                    if( (ai0ix >= ilow) && (ai0ix <= ihi) )
                    {
                        const int ivec = inc0 * (ai0ix % dim0)
                                       + dim2 * (ai1iy % dim1) 
                                              + (ai2iz % dim2);
                        assert( ivec<numpt );

                        Vector3D vc(xc,yc,zc);
                        const double r = vc.minimage(position,ll, ct.bcPoisson);
                        const double r2 = r*r;

                        const double tmp=alpha * exp(-r2*inv_rc2);
                        rhoc_[ivec] += tmp;
                        charge_ion += tmp;

                        if(r<lr)
                        {
#if 0
                            vnuc[ivec] += 
                                get_trilinval(xc,yc,zc,h0,h1,h2,position,ll,lpot);
#else
                            vnuc[ivec] += lpot.cubint(r);
#endif 
                        }
                        icount++;
                    } 

                    xc += h0;

                } // end loop over ix
           
                yc += h1;

            } // end loop over iy 

            zc += h2;

        } // end loop over iz

#if DEBUG
        if( onpe0 )
        {
            os_<<setprecision(10);
            os_<<" icount="<<icount<<endl;
            os_<<" charge compensating Ion="<<mygrid.vel()*charge_ion<<endl;
        }
#endif
        ion++;
    }

    initVcomp(ions);

    // Check compensating charges
    double  comp_rho = get_charge(rhoc_);

    if( onpe0 && ct.verbose>1 )
    {
        os_<<setprecision(8)<<fixed<<" Charge of rhoc: "<<comp_rho<<endl;
    }

#if 1
    // Rescale compensating charges
    if( onpe0 && ct.verbose>1 )
    {
        os_<<" Rescaling rhoc"<<endl;
    }
    if( ionic_charge_>0. )
    {
        double t = ionic_charge_/comp_rho;
//        my_dscal(numpt, t, rhoc_);
        MPscal(numpt, t, rhoc_);

        // Check new compensating charges
        comp_rho = get_charge(rhoc_);
    }
    if( onpe0 && ct.verbose>1 )
        os_<<" Rescaled compensating charges: "
           <<setprecision(8)<<fixed<<comp_rho<<endl;
    if(comp_rho<0.)global_exit(2);
#endif

    if( fabs(background_charge_)>0. )
    {
        double  background=background_charge_/(mygrid.gsize()*mygrid.vel());
        if( ct.bcPoisson[0]==1 && ct.bcPoisson[1]==1 && ct.bcPoisson[2]==1 )
        {
            if( onpe0 )
            {
                os_<<setprecision(12)<<scientific
                    <<"Add background charge "
                              <<background<<" to rhoc "<<endl;
            }
            for(int i=0;i<numpt;i++)
                rhoc_[i]+=background;

            // Check new compensating charges
            comp_rho = get_charge(rhoc_);
        }
    }


    electrostat_->setupRhoc(rhoc_);

    if( onpe0 && ct.verbose>3 )
        os_<<" initNuc done"<<endl;

    init_nuc_tm_.stop();
}

void MGmol::printTimers()
{  
    Control& ct = *(Control::instance());
    if( onpe0 )
    {
        os_<<setprecision(2)<<fixed<<endl;
        os_ << " Timing (real time: min_time / avg_time / max_time / min_#_calls / avg_#_calls / max_#_calls): " << endl;
        os_ << " ==================================================================== " << endl;
    }
    pb::GridFuncInterface::printTimers(os_);
    pb::GridFuncVectorInterface::printTimers(os_);
    pb::FDoperInterface::printTimers(os_);
    LocGridOrbitals::printTimers(os_);
    GridMask::printTimers(os_);
    
    sgemm_tm.print(os_);
    dgemm_tm.print(os_);
    mpgemm_tm.print(os_);
    tttgemm_tm.print(os_);
   
    ssyrk_tm.print(os_);
    dsyrk_tm.print(os_);
    mpsyrk_tm.print(os_);
    tttsyrk_tm.print(os_);
    mpdot_tm.print(os_);
    ttdot_tm.print(os_);

    dist_matrix::SubMatrices<double>::gather_tm().print(os_);
    dist_matrix::SubMatrices<double>::gather_comp_tm().print(os_);
    dist_matrix::SubMatrices<double>::gather_comm_tm().print(os_);

    dist_matrix::SparseDistMatrix<DISTMATDTYPE>::printTimers(os_);
    
    dist_matrix::DistMatrix<DISTMATDTYPE>::matgather_tm().print(os_);
    dist_matrix::DistMatrix<DISTMATDTYPE>::potrf_tm().print(os_);
    dist_matrix::DistMatrix<DISTMATDTYPE>::potri_tm().print(os_);
    
    MGmol_MPI::printTimers(os_);

    g_kbpsi_->printTimers(os_);

    get_kbpsi_tm.print(os_);
    Hamiltonian::apply_Hloc_tm().print(os_);
    computeHij_tm_.print(os_);
    Rho::printTimers(os_);
    //nonOrthoRhoKernel_tm.print(os_);
    //nonOrthoRhoKernelDiagonalBlock_tm.print(os_);
    XConGrid::get_xc_tm_.print(os_);
    get_Hpsi_and_Hij_tm_.print(os_);
    get_res_tm_.print(os_);
    comp_res_tm_.print(os_);
    vnlpsi_tm.print(os_);
    get_MLWF_tm.print(os_);
    get_NOLMO_tm.print(os_);
    Energy::eval_te_tm().print(os_);
    Electrostatic::solve_tm().print(os_);
    PoissonInterface::printTimers(os_);
    AndersonMix<LocGridOrbitals>::update_tm().print(os_);
    proj_matrices_->printTimers(os_);
    ShortSightedInverse::printTimers(os_);
    VariableSizeMatrixInterface::printTimers(os_);
    DataDistribution::printTimers(os_);
    PackedCommunicationBuffer::printTimers(os_);
    LinearSolverMatrix<lsdatatype>::printTimers(os_);
    PreconILU<pcdatatype>::printTimers(os_);
    LinearSolver::printTimers(os_);
    Table::printTimers(os_);
    LocalMatrices<MATDTYPE>::printTimers(os_);
    lrs_->printTimers(os_);
    local_cluster_->printTimers(os_);
    forces_->printTimers(os_);
    if(ct.it_algo_type == 0)
       ABPG::printTimers(os_);
    else if(ct.it_algo_type == 1)
       GrassmanLineMinimization::printTimers(os_);
    adaptLR_tm_.print(os_);
	 updateCenters_tm.print(os_);
    md_iterations_tm.print(os_);
    md_tau_tm.print(os_);
    md_moveVnuc_tm.print(os_); 
    init_nuc_tm_.print(os_); 
    md_updateMasks_tm.print(os_);        
    md_extrapolateOrbitals_tm.print(os_);         
    quench_tm.print(os_);
    evnl_tm_.print(os_);
    ions_setupInteractingIons_tm.print(os_);
    ions_setup_tm.print(os_);
    init_tm_.print(os_);
    dump_tm_.print(os_);
    setup_tm_.print(os_);
    HDFrestart::printTimers(os_);
    BlockVector<ORBDTYPE>::printTimers(os_);
    OrbitalsPreconditioning::printTimers(os_);
    MDfiles::printTimers(os_);
}

void MGmol::initKBR()
{    
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();
    Control& ct = *(Control::instance());
    Potentials& pot =hamiltonian_->potential();
    
    const double hmax=mygrid.hmax();
    vector<Species>& sp(ct.getSpecies());
    if( onpe0 && ct.verbose>0 )
        os_<<"initKBR() for "<<sp.size()<<" species, hmax="<<hmax<<endl;

    int mpi_rank;
    MPI_Comm_rank(comm_, &mpi_rank);
    int mpi_size;
    MPI_Comm_size(comm_, &mpi_size);
    
    // Distributed loop over species
    short counter=0;
    for(vector<Species>::iterator isp = sp.begin();
                                  isp!= sp.end();
                                ++isp)
    {
        if( counter%mpi_size==mpi_rank )
        {
            isp->initPotentials((bool)pot.pot_type(counter), hmax, true);
        }
        counter++;
    }
    
    //now broadcast initialized potentials to all MPI tasks
    //from task which did the work
    counter=0;
    for(vector<Species>::iterator isp = sp.begin();
                                  isp!= sp.end();
                                ++isp)
    {
        isp->syncKBP(counter%mpi_size);
        counter++;
    }
} 

double MGmol::get_evnl(const Ions& ions, LocGridOrbitals& orbitals)
{
    evnl_tm_.start();

    double val=g_kbpsi_->getEvnl(ions, orbitals, proj_matrices_);

    evnl_tm_.stop();

    return val;
}

double MGmol::getTotalEnergy(){ return total_energy_;}

void MGmol::setup()
{
    total_tm_.start();
    setup_tm_.start();

    Control& ct = *(Control::instance());
    
    if( ct.verbose>0 )
        printWithTimeStamp("MGmol::setup()...",os_);

    Mesh* mymesh = Mesh::instance();

    const int numpt=mymesh->numpt();    
    rhoc_   =new RHODTYPE[numpt];

    if( ct.verbose>0 )
        printWithTimeStamp("Setup VH...",os_);
    electrostat_=new Electrostatic(ct.lap_type,ct.bcPoisson,ct.screening_const);
    electrostat_->setup(ct.vh_init);

    rho_=new Rho();
    rho_->setVerbosityLevel(ct.verbose);
 
    int ierr=initial(); 
    if( ierr<0 ) global_exit(0); 
 
    // Write header to stdout  
    write_header();
#ifdef USE_MPI
    if( ct.verbose>5 )
    {
        Mesh* mymesh = Mesh::instance();
        const pb::PEenv& myPEenv=mymesh->peenv();
        myPEenv.printPEnames(os_);
    }
#endif
   
    if( ct.verbose>0 )
        printWithTimeStamp("MGmol::setup done...",os_);
     
    setup_tm_.stop();
}

void MGmol::cleanup()
{
    closing_tm_.start();

    Mesh* mymesh = Mesh::instance();
    const pb::PEenv& myPEenv=mymesh->peenv();
    Control& ct = *(Control::instance());

    printTimers();

    // Save data to restart file  
    if( ct.out_restart_info>0 
     && ct.atoms_dyn==0 )
    {
        const pb::Grid& mygrid  = mymesh->grid();
        unsigned gdim[3]={mygrid.gdim(0),mygrid.gdim(1),mygrid.gdim(2)};
    
        // create restart file
        string filename(string(ct.out_restart_file));
        filename+="0";
        HDFrestart h5restartfile(filename, myPEenv, gdim, ct.out_restart_file_type);

        int ierr=write_hdf5(h5restartfile, rho_->rho_, *ions_, *current_orbitals_, *lrs_);
        
        if( ierr<0 )
            os_<<"WARNING: writing restart data failed!!!"<<endl;
    }

    MPI_Barrier( comm_ );
    closing_tm_.stop();
    total_tm_.stop();
    if( onpe0 ){
        os_ << " ******************************************************************** "  << endl;        
    }
        closing_tm_.print(os_);
        total_tm_.print(os_);
}

void MGmol::projectOutKernel(LocGridOrbitals& phi)
{
    assert( aomm_!=0 );
    aomm_->projectOut(phi);
}

void MGmol::setGamma(const pb::Lap<ORBDTYPE>& lapOper, const Potentials& pot)
{
    assert( orbitals_precond_!=0 );
    
    Control& ct = *(Control::instance());

    orbitals_precond_->setGamma(lapOper,pot,ct.getMGlevels(),proj_matrices_);
}

void MGmol::precond_mg(LocGridOrbitals& phi)
{
    assert( orbitals_precond_!=0 );
    
    orbitals_precond_->precond_mg(phi);
}

double MGmol::computeResidual(LocGridOrbitals& orbitals,
            LocGridOrbitals& work_orbitals,
            LocGridOrbitals& res,
            const bool print_residual,
            const bool norm_res)

{
    assert( orbitals.getIterativeIndex()>=0 );
    
    comp_res_tm_.start();
    //os_<<"computeResidual()"<<endl;

    Control& ct ( *(Control::instance()) );

    proj_matrices_->computeInvB();

    Potentials& pot =hamiltonian_->potential();
    pb::Lap<ORBDTYPE>* lapop=hamiltonian_->lapOper();

    setGamma(*lapop,pot);

    // get H*psi stored in work_orbitals.psi
    // and psi^T H psi in Hij
    getHpsiAndTheta(*ions_, orbitals, work_orbitals);
    
    double norm2Res=
    computeConstraintResidual(orbitals, 
                                      work_orbitals,
                                      res, 
                                      print_residual,
                                      norm_res);


    if( ct.isSpreadFunctionalEnergy() )
        addResidualSpreadPenalty(orbitals,res);

    comp_res_tm_.stop();
    
    return norm2Res;
}

//////////////////////////////////////////////////////////////////////////////
// compute res using psi and hpsi
void MGmol::computeResidualUsingHPhi(LocGridOrbitals& psi,
                        const LocGridOrbitals& hphi,
                        LocGridOrbitals& res,
                        const bool applyB)
{
    assert( psi.isCompatibleWith(hphi) );
    assert( psi.isCompatibleWith(res) );
    
    get_res_tm_.start();

    proj_matrices_->updateSubMatT();
    
    SquareLocalMatrices<MATDTYPE>& localT( proj_matrices_->getLocalT() );

    pb::Lap<ORBDTYPE>* lapop=hamiltonian_->lapOper();
    const int ncolors=psi.chromatic_number();
    if( ncolors>0 )
    {
        // compute B*psi and store in tmp
        ORBDTYPE* old_storage=0;
        vector<ORBDTYPE> tmp;
        if( applyB )
        {
            const int ld=psi.getLda();
            tmp.resize(ld*ncolors);
            psi.setDataWithGhosts();
            psi.trade_boundaries();
            for(int i=0;i<ncolors;i++)
            {
                lapop->rhs(psi.getFuncWithGhosts(i),&tmp[0]+ld*i);
            }

            // psi points to tmp temporarily
            old_storage=psi.getPsi(0);    
            psi.set_storage(&tmp[0]);
        }
        
        // get B*phi*theta and store it in res in [Ry]
        // (even if Ritz functions mode)
        psi.multiplyByMatrix(0, ncolors, localT, res);
        
        if( applyB )
        {
            psi.set_storage(old_storage);
        }
        
        // res = (B*phi*theta - H*phi) in [Ry]
        res.axpy(-1.,hphi);
    }
 
#if 0
    dist_matrix::DistMatrix<DISTMATDTYPE> RPsi(psi.product(res));
    if( onpe0 )
        os_<<" matrix R**T * Psi\n";
    RPsi.print(os_,0,0,5,5);
    if(onpe0)os_<<"trace = "<<RPsi.trace()<<endl;;    
#endif

    get_res_tm_.stop();
}

//////////////////////////////////////////////////////////////////////////////

double MGmol::computeConstraintResidual(LocGridOrbitals& orbitals, 
                         const LocGridOrbitals& hphi, 
                         LocGridOrbitals& res, 
                         const bool print_residual,
                         const bool compute_norm_res)
{
    Control& ct ( *(Control::instance()) );

    computeResidualUsingHPhi(orbitals, hphi, res, ct.Mehrstellen());

    res.applyCorrMask();

    double normRes=-1.;
    if( print_residual || compute_norm_res )
    {
        if( ct.checkMaxResidual() )
        {
            normRes=0.5*res.maxAbsValue(); // factor 0.5 for Rydberg to Hartree conversion
            if(onpe0 && print_residual)
                os_<<setprecision(2)<<scientific
                    <<"max. |Residual_ij| = "<<normRes<<endl; 
        }
        else
        {
            double norm2Res=res.dotProduct(res);
            normRes=0.5*sqrt(norm2Res); // factor 0.5 for Rydberg to Hartree conversion
            if(onpe0 && print_residual)
                os_<<setprecision(2)<<scientific
//                   <<"|| Relative residual || = "<<normRes/ct.getNel()<<endl; 
                   <<"|| Residual || = "<<normRes<<endl; 
        }
    }

    return normRes;
}

//////////////////////////////////////////////////////////////////////////////
// Get preconditioned residual in res_orbitals
double MGmol::computePrecondResidual(LocGridOrbitals& phi,
                        LocGridOrbitals& hphi,
                        LocGridOrbitals& res,
                        Ions& ions,
                        KBPsiMatrixInterface* kbpsi, 
                        const bool print_residual,
                        const bool norm_res)

{
    Control& ct = *(Control::instance());

    proj_matrices_->computeInvB();

    Potentials& pot =hamiltonian_->potential();
    pb::Lap<ORBDTYPE>* lapop=hamiltonian_->lapOper();

    setGamma(*lapop,pot);

    // get H*psi stored in hphi
    // and psi^T H psi in Hij
    getHpsiAndTheta(ions, phi, hphi, kbpsi);

    double norm2Res=
    computeConstraintResidual(phi, 
                          hphi,
                          res,
                          print_residual,
                          norm_res);

    if( (ct.getPrecondType()%10)==0 && ct.getMGlevels()>=0 )
    {
        // PRECONDITIONING 
        // compute the preconditioned steepest descent direction
        // -> res
        orbitals_precond_->precond_mg(res);
    }

    //if( ct.isSpreadFunctionalActive() )addResidualSpreadPenalty(phi,res);
    
    return norm2Res;
} 

//    Function to update potentials vh and vxc:
//    
//    The new potentials are computed as a linear combination 
//    of the old ones (input "vh" and "vxc") and the ones 
//    corresponding to the input "rho".
void MGmol::update_pot(const pb::GridFunc<POTDTYPE>& vh_init,
                const Ions& ions)
{
    electrostat_->setupInitialVh(vh_init);
    update_pot(ions);
}

void MGmol::update_pot(const Ions& ions)
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )os_<<"Update potentials"<<endl;
#endif

    Control& ct = *(Control::instance());
    Potentials& pot =hamiltonian_->potential();
    
    // Update exchange-correlation potential
    xcongrid_->update();

    // Generate new hartree potential
    electrostat_->computeVh(ions, *rho_, pot);

    const bool flag_mixing = ( fabs(ct.mix_pot-1.)>1.e-3 );

    // evaluate potential correction
    if( flag_mixing )
    {
        pot.delta_v(rho_->rho_);
        pot.update(ct.mix_pot);
    }
    else
        pot.update(rho_->rho_);

#if 0
    if ( onpe0 )
    {
        os_<<setprecision(3);
        os_<<" <rho dv>/nb_ions = "<<pot.scf_dvrho()/ions.num_ions()<<endl;
    }
#endif
}


void MGmol::addResidualSpreadPenalty(LocGridOrbitals& phi,
                                     LocGridOrbitals& res)
{
    assert( spread_penalty_!=0 );
    
    spread_penalty_->addResidual(phi,res);
}

