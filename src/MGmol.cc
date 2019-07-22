// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

#include "global.h"

#include "ABPG.h"
#include "AOMMprojector.h"
#include "AndersonMix.h"
#include "ConstraintSet.h"
#include "Control.h"
#include "DFTsolver.h"
#include "DMStrategyFactory.h"
#include "DistMatrix.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "EnergySpreadPenalty.h"
#include "FDoper.h"
#include "FIRE.h"
#include "Forces.h"
#include "GrassmanLineMinimization.h"
#include "GridFunc.h"
#include "HDFrestart.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "KBPsiMatrixSparse.h"
#include "LBFGS.h"
#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "MDfiles.h"
#include "MGmol.h"
#include "MLWFTransform.h"
#include "MPIdata.h"
#include "MasksSet.h"
#include "Mesh.h"
#include "OrbitalsPreconditioning.h"
#include "PackedCommunicationBuffer.h"
#include "PoissonInterface.h"
#include "Potentials.h"
#include "Power.h"
#include "Preconditioning.h"
#include "ProjectedMatricesMehrstellen.h"
#include "ProjectedMatricesSparse.h"
#include "Rho.h"
#include "SparseDistMatrix.h"
#include "SpreadPenalty.h"
#include "SpreadPenaltyVolume.h"
#include "SpreadsAndCenters.h"
#include "SubMatrices.h"
#include "SubspaceProjector.h"
#include "XConGrid.h"
#include "XCfunctionalFactory.h"
#include "DistMatrix2SquareLocalMatrices.h"
#include "SP2.h"
#include "manage_memory.h"

namespace mgmol
{
std::ostream* out = NULL;
}

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

#include "Signal.h"
set<int> Signal::recv_;

template <class T>
MGmol<T>::MGmol(MPI_Comm comm, std::ostream& os) : os_(os)
{
    comm_ = comm;

    constraints_ = new ConstraintSet();

    mgmol::out = &os;

    geom_optimizer_ = 0;
    lrs_            = 0;
    local_cluster_  = 0;
    proj_matrices_  = 0;
    dm_strategy_    = 0;

    h5f_file_ = 0;

    aomm_ = 0;

    spreadf_ = 0;

    spread_penalty_ = 0;

    data_distributor_ = new BasicDataDistributors();

    orbitals_precond_ = 0;

    forces_ = 0;

    energy_ = 0;
}

template <class T>
MGmol<T>::~MGmol()
{
    Control& ct = *(Control::instance());
    delete electrostat_;
    delete rho_;
    delete constraints_;
    if (!ct.short_sighted)
    {
        delete remote_tasks_DistMatrix_;
        remote_tasks_DistMatrix_ = 0;
    }
    delete xcongrid_;
    delete energy_;
    if (hamiltonian_ != 0) delete hamiltonian_;
    if (geom_optimizer_ != 0) delete geom_optimizer_;

    delete currentMasks_;
    delete corrMasks_;

    if (aomm_ != 0) delete aomm_;

    delete current_orbitals_;
    delete ions_;
    delete g_kbpsi_;

    delete proj_matrices_;
    delete lrs_;
    if (local_cluster_ != 0) delete local_cluster_;

    if (h5f_file_ != 0) delete h5f_file_;

    if (spreadf_ != 0) delete spreadf_;

    if (spread_penalty_ != 0) delete spread_penalty_;

    delete data_distributor_;

    delete forces_;
    if (dm_strategy_ != 0) delete dm_strategy_;
}

template <class T>
int MGmol<T>::initial()
{
    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Mesh* mymesh    = Mesh::instance();

    assert(ct.numst >= 0);

    if (ct.verbose > 0) printWithTimeStamp("MGmol<T>::initial()...", os_);

    init_tm_.start();

    const pb::Grid& mygrid = mymesh->grid();

    hamiltonian_->setup(mygrid, ct.lap_type);

    pb::Lap<ORBDTYPE>* lapop = ct.Mehrstellen() ? hamiltonian_->lapOper() : 0;
    g_kbpsi_                 = new KBPsiMatrixSparse(lapop);

    check_anisotropy();

    if (ct.verbose > 0)
        printWithTimeStamp(
            "MGmol<T>::initial(), create ProjectedMatrices...", os_);

    // If not an initial run read data from files
    if (ct.restart_info > 2 && ct.isLocMode())
    {
        string name = "ExtrapolatedFunction";
        if (ct.verbose > 0)
            printWithTimeStamp(
                "read LRs from ExtrapolatedFunction database...", os_);
        int n = read_restart_lrs(*h5f_file_, name);
        if (n == 0)
        {
            if (ct.verbose > 0)
                printWithTimeStamp("read LRs from Function database...", os_);
            name = "Function";
            n    = read_restart_lrs(*h5f_file_, name);
        }
        if (n < 0) return n;

        if (n > 0) lrs_->setup();
    }

    double dlrsmin
        = lrs_->computeMinDistBetweenLocalPairs(cout, (ct.verbose > 2));
    if (dlrsmin < 1.e-3)
    {
        cout << "WARNING: Min. distance between LR centers is " << dlrsmin
             << "!!!" << endl;
    }

    // initialize and setup load balancing object
    // set number of iterations to 10.
    if (ct.load_balancing_alpha > 0.0)
    {
        local_cluster_ = new ClusterOrbitals(lrs_);
        local_cluster_->setup();
        local_cluster_->computeClusters(ct.load_balancing_max_iterations);
    }
    // initialize data distribution objects
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    data_distributor_->initialize(lrs_, myPEenv, domain);

    bool with_spin = (mmpi.nspin() > 1);
    if (ct.Mehrstellen())
        proj_matrices_ = new ProjectedMatricesMehrstellen(ct.numst, with_spin);
    else if (ct.short_sighted)
        proj_matrices_
            = new ProjectedMatricesSparse(ct.numst, lrs_, local_cluster_);
    else
        proj_matrices_ = new ProjectedMatrices(ct.numst, with_spin);

    forces_ = new Forces<T>(hamiltonian_, rho_, proj_matrices_);

    if (ct.verbose > 0)
        printWithTimeStamp("MGmol<T>::initial(), create MasksSet...", os_);

    currentMasks_ = new MasksSet(false, ct.getMGlevels());
    currentMasks_->setup(*lrs_);

    corrMasks_ = new MasksSet(true, 0);
    corrMasks_->setup(*lrs_);

    BlockVector<ORBDTYPE>::setOverAllocateFactor(
        ct.orbitalsOverallocateFactor());

    if (ct.verbose > 0)
        printWithTimeStamp("MGmol<T>::initial(), create T...", os_);

    if (!ct.short_sighted)
    {
        printWithTimeStamp(
            "MGmol<T>::initial(), create MatricesBlacsContext and misc...", os_);

        dist_matrix::DistMatrix<DISTMATDTYPE> tmp("tmp", ct.numst, ct.numst);
        remote_tasks_DistMatrix_
            = new dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>(tmp);
        ProjectedMatrices::registerRemoteTasksDistMatrix(
            remote_tasks_DistMatrix_);
        remote_tasks_DistMatrix_ptr_ = remote_tasks_DistMatrix_;
    }

    current_orbitals_ = new T("Primary", mygrid,
        mymesh->subdivx(), ct.numst,
        ct.bc, proj_matrices_, lrs_, currentMasks_, corrMasks_, local_cluster_,
        true);

    increaseMemorySlotsForOrbitals();

    Potentials& pot            = hamiltonian_->potential();
    pb::Lap<ORBDTYPE>* lapOper = hamiltonian_->lapOper();

    // If not an initial run read data from files
    if (ct.restart_info > 0)
    {
        md_time_      = h5f_file_->getMDTimeFromFile();
        md_iteration_ = h5f_file_->getMDstepFromFile();

        int ierr = read_restart_data(*h5f_file_, *rho_, *current_orbitals_);
        if (ierr < 0) return -1;

        if (ct.restart_info > 2)
        {
            current_orbitals_->applyMask();
        }
    }
    else
    {
        md_time_      = 0.;
        md_iteration_ = 1;
    }

    mmpi.barrier();
    if (ct.verbose > 0) current_orbitals_->printChromaticNumber(os_);

    pot.initBackground(*ions_);

    // Random initialization of the wavefunctions
    if (ct.restart_info <= 2)
    {
        if (ct.verbose > 0)
            printWithTimeStamp("MGmol<T>::initial(), init wf and masks...", os_);

        // Make temp mask for initial random wave functions
        if (ct.init_loc == 1)
        {
            float cut_init = ct.initRadius();
            assert(cut_init > 0.);
            currentMasks_->initialize(*lrs_, 0, cut_init);
        }

        // Initialize states
        current_orbitals_->initWF(*lrs_);

        // initialize masks again
        if (ct.init_loc == 1)
        {
            currentMasks_->initialize(*lrs_, 0);
        }
    }

    // Initialize the radial potential stuff
    if (ct.verbose > 0) printWithTimeStamp("initKBR()...", os_);
    initKBR();

    ions_->setup();

    double d = ions_->computeMinLocalSpacing();
    if (d < 1.e-3)
    {
        cerr << "ERROR: min. distance between ions is smaller than 1.e-3!!!\n";
    }

    // Initialize the nuclear local potential and the compensating charges
    if (ct.verbose > 0) printWithTimeStamp("initNuc()...", os_);
    initNuc(*ions_);

    // initialize Rho
    if (ct.verbose > 0) printWithTimeStamp("Initialize Rho...", os_);
    if (ct.restart_info <= 1)
        proj_matrices_->setDMuniform(
            ct.getNelSpin(), current_orbitals_->getIterativeIndex());

    rho_->setup(ct.getOrbitalsType(), current_orbitals_->getOverlappingGids());

    if (ct.restart_info <= 1)
    {
        if (ct.init_type == 0) // random functions
        {
            rho_->init(pot.rho_comp());
        }
        else
        {
            rho_->update(*current_orbitals_);
        }
    }

    if (ct.verbose > 0) printWithTimeStamp("Initialize XC functional...", os_);
    xcongrid_ = XCfunctionalFactory<T>::create(
                    ct.xctype, mmpi.nspin(), *rho_, pot);
    assert( xcongrid_ != 0 );

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

    if (ct.numst > 0)
    {
        if (ct.verbose > 0)
            printWithTimeStamp("Compute initial Gram matrix...", os_);
        current_orbitals_->computeGramAndInvS(ct.verbose);
        if (ct.verbose > 0)
            printWithTimeStamp("Compute initial B matrix...", os_);
        current_orbitals_->computeBAndInvB(*lapOper);
        if (ct.verbose > 0)
            printWithTimeStamp("Compute initial condition number...", os_);
        current_orbitals_->checkCond(100000., ct.AtomsMove());
    }

    if (ct.verbose > 0) printWithTimeStamp("Setup kbpsi...", os_);
    g_kbpsi_->setup(*ions_, *current_orbitals_);

    if (ct.restart_info == 0)
    {
        if (ct.verbose > 0) printWithTimeStamp("update_pot...", os_);
        update_pot(*ions_);
    }

    if (ct.isLocMode() || ct.isSpreadFunctionalActive())
    {
        Vector3D origin(mygrid.origin(0), mygrid.origin(1), mygrid.origin(2));
        Vector3D ll(mygrid.ll(0), mygrid.ll(1), mygrid.ll(2));
        spreadf_ = new SpreadsAndCenters<T>(origin, ll);
    }

    bool energy_with_spread_penalty = false;
    if (ct.isSpreadFunctionalActive())
    {
        if (ct.isSpreadFunctionalVolume())
        {
            spread_penalty_ = new SpreadPenaltyVolume<T>(spreadf_,
                ct.spreadPenaltyTarget(), ct.spreadPenaltyAlphaFactor(),
                ct.spreadPenaltyDampingFactor());
        }
        else if (ct.isSpreadFunctionalEnergy())
        {
            energy_with_spread_penalty = true;
            spread_penalty_            =
                new EnergySpreadPenalty<T>(spreadf_,
                    ct.spreadPenaltyTarget(), ct.spreadPenaltyAlphaFactor());
        }
        else
            spread_penalty_ = new SpreadPenalty<T>(spreadf_,
                ct.spreadPenaltyTarget(), ct.spreadPenaltyAlphaFactor(),
                ct.spreadPenaltyDampingFactor());
    }

    SpreadPenaltyInterface<T>* spread_penalty
        = energy_with_spread_penalty ? spread_penalty_ : 0;
    energy_ = new Energy<T>(
        mygrid, *ions_, pot, *electrostat_, *rho_, *xcongrid_, spread_penalty);

    if (ct.verbose > 0) printWithTimeStamp("Setup matrices...", os_);

    updateHmatrix(*current_orbitals_, *ions_);

    // HMVP algorithm requires that H is initialized
    dm_strategy_ = DMStrategyFactory<T>::create(
        comm_, os_, *ions_, rho_, energy_,
        electrostat_, this, proj_matrices_, current_orbitals_);

    // theta = invB * Hij
    proj_matrices_->updateThetaAndHB();

    dm_strategy_->initialize();

    if (ct.verbose > 1)
    {
        proj_matrices_->printMatrices(os_);
        proj_matrices_->printDM(os_);
    }

    init_tm_.stop();

    if (ct.verbose > 0) printWithTimeStamp("Initialization done...", os_);

    return 0;
} // initial()

template <class T>
void MGmol<T>::run()
{
    total_tm_.start();

    setup();

    Control& ct = *(Control::instance());

    double eks = 0.;

    if (ct.verbose > 0) printWithTimeStamp("Run...", os_);
    // Dispatch to the method chosen
    switch (ct.AtomsDynamic())
    {
        case AtomsDynamicType::Quench:
            quench(current_orbitals_, *ions_, ct.max_electronic_steps, 20, eks);

            // Forces for the last states
            force(*current_orbitals_, *ions_);

            constraints_->addConstraints(*ions_);

            constraints_->setup(*ions_);

            constraints_->printConstraintsForces(os_);

            constraints_->projectOutForces(20);

            if ((ions_->getNumIons() <= 1024 || ct.verbose > 2)
                && ct.verbose > 0)
                ions_->printForcesGlobal(os_);

            finalEnergy();

            printMM();

            break;

        case AtomsDynamicType::MD:
            md(&current_orbitals_, *ions_);
            // finalEnergy();
            break;

        case AtomsDynamicType::LBFGS:
            lbfgsrlx(&current_orbitals_, *ions_);
            // finalEnergy();
            break;

        case AtomsDynamicType::FIRE:
            runfire(&current_orbitals_, *ions_);
            // finalEnergy();
            break;

        default:
            (*MPIdata::serr) << "run: Undefined MD method" << endl;
    }

    cleanup();
} // run()

template <class T>
void MGmol<T>::finalEnergy()
{
    // Get the total energy
    const double ts = 0.5 * proj_matrices_->computeEntropy(); // in [Ha]
    total_energy_   = energy_->evaluateTotal(
        ts, proj_matrices_, *current_orbitals_, 2, os_);
}

template <class T>
void MGmol<T>::printMM()
{
    Control& ct = *(Control::instance());
    if (ct.tmatrices == 1)
    {
        ofstream tfile("s.mm", ios::out);
        proj_matrices_->printGramMM(tfile);
        ofstream tfileh("h.mm", ios::out);
        ProjectedMatrices* projmatrices
            = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
        assert(projmatrices != 0);
        projmatrices->printHamiltonianMM(tfileh);
    }
}

/* Writes out header information */
template <class T>
void MGmol<T>::write_header()
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    time_t tt;

    time(&tt);
    char* timeptr = ctime(&tt);

    Control& ct = *(Control::instance());
    Potentials& pot = hamiltonian_->potential();

    if (onpe0)
    {

        os_ << "//////////////////////////////////////////////////////////"
            << endl;
        os_ << endl;
#ifdef GITHASH
#define xstr(x) #x
#define LOGGIT(x) os_ << " MGmol: git_hash " << xstr(x) << endl;
        LOGGIT(GITHASH);
        os_ << endl;
#endif
        os_ << " Compiled: " << __DATE__ << ", " << __TIME__ << endl;
        os_ << " Real-space finite difference ab initio calculations\n";
        os_ << endl;
        os_ << " authors: J.-L. Fattebert, Oak Ridge National Laboratory\n";
        os_ << "          D. Osei-Kuffuor, Lawrence Livermore National "
               "Laboratory\n";
        os_ << "          I.S. Dunn, Columbia University\n\n";
        os_ << "//////////////////////////////////////////////////////////"
            << endl;

        os_ << endl << endl << description << endl;
        os_ << " Run started at " << timeptr << endl;

        os_ << " Orbitals precision (in bytes): " << sizeof(ORBDTYPE) << endl
            << endl;

        Potentials& pot = hamiltonian_->potential();
        pot.writeNames(os_);

        mymesh->print(os_);

        pb::Lap<ORBDTYPE>* lapop = hamiltonian_->lapOper();
        os_ << " Laplacian Discretization: " << lapop->name() << endl;

#ifdef SMP_NODE
        os_ << " " << omp_get_max_threads() << " thread"
            << (omp_get_max_threads() > 1 ? "s " : " ");
        os_ << "active" << endl << endl;
#endif

        os_ << " ScaLapack block size: "
            << dist_matrix::DistMatrix<DISTMATDTYPE>::getBlockSize() << endl;

        if (!ct.short_sighted)
        {
            MatricesBlacsContext& mbc(MatricesBlacsContext::instance());
            mbc.print(os_);
        }
    } // onpe0

    if (current_orbitals_ != NULL && ct.verbose > 0)
    {
        current_orbitals_->printNumst(os_);
        current_orbitals_->printChromaticNumber(os_);
    }

    int nions = ions_->getNumIons();
    if (onpe0)
    {
        os_ << " Number of ions     = " << nions << endl;
        os_ << " Total charge in cell = " << pot.getChargeInCell() << endl;

        ct.print(os_);

        os_ << endl << endl << " Atomic species information:" << endl << endl;
        os_ << fixed << setprecision(5);

        const std::vector<Species>& sp(ct.getSpecies());
        for (int idx = 0; idx < (int)sp.size(); idx++)
        {
            const Species& sp_ion(sp[idx]);

            os_ << " Species #" << idx + 1 << ": ";
            sp_ion.print(os_);

            os_ << " dim_l        = " << sp_ion.dim_l()
                << "    ->Diameter    = " << sp_ion.dim_l() * mygrid.hmin()
                << "[bohr]" << endl;
            os_ << " dim_nl      = " << sp_ion.dim_nl()
                << "    ->Diameter    = " << sp_ion.dim_nl() * mygrid.hmin()
                << "[bohr]" << endl;
        }
        if (ct.short_sighted)
        {
            os_ << endl;
            os_ << " Short_Sighted Solver Parameters: " << endl << endl;
            os_ << " Accelerator = Flexible GMRES(" << ct.fgmres_kim << ")"
                << endl;
            if (ct.ilu_type == 0)
                os_ << " Preconditioner = ilu" << ct.ilu_lof << endl;
            else if (ct.ilu_type == 1)
                os_ << " Preconditioner = modified ilu" << endl;
            else if (ct.ilu_type == 2)
                os_ << " Preconditioner = standard ilu" << endl;
            else
                os_ << " Preconditioner = ilu0 preconditioner" << endl;
            os_ << scientific
                << " fgmres convergence tolerance = " << ct.fgmres_tol << endl;
            os_ << " fgmres max. iterations = " << ct.fgmres_maxits << endl;
            os_ << scientific
                << " drop tolerance for (standard/ modified) ilu = "
                << ct.ilu_droptol << endl;
            os_ << " max. fill-in for (standard/ modified) ilu = "
                << ct.ilu_maxfil << endl;
        }
    } // onpe0

    // Write out the ionic postions and displacements
    ions_->printPositions(os_);

    if (current_orbitals_ != NULL && ct.verbose > 3) lrs_->printAllRegions(os_);
}

template <class T>
void MGmol<T>::global_exit(int i)
{
#ifdef USE_MPI
    MPI_Abort(comm_, i);
#else
    abort();
#endif
}

template <class T>
void MGmol<T>::check_anisotropy()
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    if (mygrid.anisotropy() > 1.2)
    {
        write_header();
        if (onpe0)
            (*MPIdata::serr) << " hmax=" << mygrid.hmax()
                             << ", hmin=" << mygrid.hmin() << endl;
        (*MPIdata::serr) << "init: Anisotropy too large: "
                         << mygrid.anisotropy() << endl;
        global_exit(2);
    }
}

template <class T>
void MGmol<T>::printEigAndOcc()
{
    Control& ct = *(Control::instance());
    if (!(ct.fullyOccupied()
            && ct.getOrbitalsType() != OrbitalsType::Eigenfunctions
            && ct.occupationWidthIsZero())
        && onpe0)
    {
        ProjectedMatrices* projmatrices
            = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
        assert(projmatrices);

        projmatrices->printEigenvalues(os_);
        projmatrices->printOccupations(os_);
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

template <class T>
void MGmol<T>::initNuc(Ions& ions)
{
    init_nuc_tm_.start();

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 1) os_ << "Init vnuc and rhoc..." << endl;

    Potentials& pot = hamiltonian_->potential();

    pot.initialize(ions);

    // Check compensating charges
    double comp_rho = pot.getCharge(pot.rho_comp());

    if (onpe0 && ct.verbose > 1)
    {
        os_ << setprecision(8) << fixed << " Charge of rhoc: " << comp_rho
            << endl;
    }

#if 1
    pot.rescaleRhoComp();
#endif

    pot.addBackgroundToRhoComp();

    electrostat_->setupRhoc(pot.rho_comp());

    if (onpe0 && ct.verbose > 3) os_ << " initNuc done" << endl;

    init_nuc_tm_.stop();
}

template <class T>
void MGmol<T>::printTimers()
{
    Control& ct = *(Control::instance());
    if (onpe0)
    {
        os_ << setprecision(2) << fixed << endl;
        os_ << " Timing (real time: min_time / avg_time / max_time / "
               "min_#_calls / avg_#_calls / max_#_calls): "
            << endl;
        os_ << " =============================================================="
               "====== "
            << endl;
    }
    pb::GridFuncInterface::printTimers(os_);
    pb::GridFuncVectorInterface::printTimers(os_);
    pb::FDoperInterface::printTimers(os_);
    T::printTimers(os_);
    SinCosOps<T>::printTimers(os_);
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

    dist_matrix::SubMatrices<double>::printTimers(os_);

    DistMatrix2SquareLocalMatrices::printTimers(os_);

    dist_matrix::SparseDistMatrix<DISTMATDTYPE>::printTimers(os_);

    dist_matrix::DistMatrix<DISTMATDTYPE>::printTimers(os_);

    MGmol_MPI::printTimers(os_);

    g_kbpsi_->printTimers(os_);

    get_kbpsi_tm.print(os_);
    Hamiltonian<T>::apply_Hloc_tm().print(os_);
    computeHij_tm_.print(os_);
    Rho<T>::printTimers(os_);
    XConGrid::get_xc_tm_.print(os_);
    get_Hpsi_and_Hij_tm_.print(os_);
    get_res_tm_.print(os_);
    comp_res_tm_.print(os_);
    vnlpsi_tm.print(os_);
    get_MLWF_tm.print(os_);
    get_NOLMO_tm.print(os_);
    Energy<T>::eval_te_tm().print(os_);
    Electrostatic::solve_tm().print(os_);
    PoissonInterface::printTimers(os_);
    AndersonMix<T>::update_tm().print(os_);
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
    Power<LocalVector<double>, SquareLocalMatrices<double>>::printTimers(os_);
    SP2::printTimers(os_);
    lrs_->printTimers(os_);
    local_cluster_->printTimers(os_);
    forces_->printTimers(os_);
    if (ct.OuterSolver() == OuterSolverType::ABPG)
        ABPG<T>::printTimers(os_);
    else if (ct.OuterSolver() == OuterSolverType::NLCG)
        GrassmanLineMinimization<T>::printTimers(os_);
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
    OrbitalsPreconditioning<T>::printTimers(os_);
    MDfiles::printTimers(os_);
}

template <class T>
void MGmol<T>::initKBR()
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    Control& ct            = *(Control::instance());
    Potentials& pot        = hamiltonian_->potential();

    const double hmax = mygrid.hmax();
    vector<Species>& sp(ct.getSpecies());
    if (onpe0 && ct.verbose > 0)
        os_ << "initKBR() for " << sp.size() << " species, hmax=" << hmax
            << endl;

    int mpi_rank;
    MPI_Comm_rank(comm_, &mpi_rank);
    int mpi_size;
    MPI_Comm_size(comm_, &mpi_size);

    // Distributed loop over species
    short counter = 0;
    for (vector<Species>::iterator isp = sp.begin(); isp != sp.end(); ++isp)
    {
        if (counter % mpi_size == mpi_rank)
        {
            isp->initPotentials((bool)pot.pot_type(counter), hmax, true);
        }
        counter++;
    }

    // now broadcast initialized potentials to all MPI tasks
    // from task which did the work
    counter = 0;
    for (vector<Species>::iterator isp = sp.begin(); isp != sp.end(); ++isp)
    {
        isp->syncKBP(counter % mpi_size);
        counter++;
    }
}

template <class T>
double MGmol<T>::get_evnl(const Ions& ions, T& orbitals)
{
    evnl_tm_.start();
    Control& ct = *(Control::instance());

    double val;
    if (ct.short_sighted)
    {
        ProjectedMatricesSparse* projmatrices
            = dynamic_cast<ProjectedMatricesSparse*>(proj_matrices_);
        assert(projmatrices);

        val = g_kbpsi_->getEvnl(ions, orbitals, projmatrices);
    }
    else
    {
        ProjectedMatrices* projmatrices
            = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
        assert(projmatrices);

        val = g_kbpsi_->getEvnl(ions, orbitals, projmatrices);
    }

    evnl_tm_.stop();

    return val;
}

template <class T>
double MGmol<T>::getTotalEnergy() { return total_energy_; }

template <class T>
void MGmol<T>::setup()
{
    total_tm_.start();
    setup_tm_.start();

    Control& ct = *(Control::instance());

    if (ct.verbose > 0) printWithTimeStamp("MGmol<T>::setup()...", os_);

    if (ct.verbose > 0) printWithTimeStamp("Setup VH...", os_);
    electrostat_
        = new Electrostatic(ct.lap_type, ct.bcPoisson, ct.screening_const);
    electrostat_->setup(ct.vh_init);

    rho_ = new Rho<T>();
    rho_->setVerbosityLevel(ct.verbose);

    int ierr = initial();
    if (ierr < 0) global_exit(0);

    // Write header to stdout
    write_header();
#ifdef USE_MPI
    if (ct.verbose > 5)
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        myPEenv.printPEnames(os_);
    }
#endif

    if (ct.verbose > 0) printWithTimeStamp("MGmol<T>::setup done...", os_);

    setup_tm_.stop();
}

template <class T>
void MGmol<T>::cleanup()
{
    closing_tm_.start();

    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    Control& ct              = *(Control::instance());

    printTimers();

    // Save data to restart file
    if (ct.out_restart_info > 0 && !ct.AtomsMove())
    {
        const pb::Grid& mygrid = mymesh->grid();
        unsigned gdim[3] = { mygrid.gdim(0), mygrid.gdim(1), mygrid.gdim(2) };

        // create restart file
        string filename(string(ct.out_restart_file));
        filename += "0";
        HDFrestart h5restartfile(
            filename, myPEenv, gdim, ct.out_restart_file_type);

        int ierr = write_hdf5(
            h5restartfile, rho_->rho_, *ions_, *current_orbitals_, *lrs_);

        if (ierr < 0) os_ << "WARNING: writing restart data failed!!!" << endl;
    }

    MPI_Barrier(comm_);
    closing_tm_.stop();
    total_tm_.stop();
    if (onpe0)
    {
        os_ << " **************************************************************"
               "****** "
            << endl;
    }
    closing_tm_.print(os_);
    total_tm_.print(os_);
}

template <>
void MGmol<LocGridOrbitals>::projectOutKernel(LocGridOrbitals& phi)
{
    assert(aomm_ != 0);
    aomm_->projectOut(phi);
}

template <class T>
void MGmol<T>::projectOutKernel(T& phi)
{
    (void)phi;

    return;
}

template <class T>
void MGmol<T>::setGamma(const pb::Lap<ORBDTYPE>& lapOper, const Potentials& pot)
{
    assert(orbitals_precond_ != 0);

    Control& ct = *(Control::instance());

    orbitals_precond_->setGamma(lapOper, pot, ct.getMGlevels(), proj_matrices_);
}

template <class T>
void MGmol<T>::precond_mg(T& phi)
{
    assert(orbitals_precond_ != 0);

    orbitals_precond_->precond_mg(phi);
}

template <class T>
double MGmol<T>::computeResidual(T& orbitals,
    T& work_orbitals, T& res,
    const bool print_residual, const bool norm_res)

{
    assert(orbitals.getIterativeIndex() >= 0);

    comp_res_tm_.start();
    // os_<<"computeResidual()"<<endl;

    Control& ct(*(Control::instance()));

    proj_matrices_->computeInvB();

    Potentials& pot          = hamiltonian_->potential();
    pb::Lap<ORBDTYPE>* lapop = hamiltonian_->lapOper();

    setGamma(*lapop, pot);

    // get H*psi stored in work_orbitals.psi
    // and psi^T H psi in Hij
    getHpsiAndTheta(*ions_, orbitals, work_orbitals);

    double norm2Res = computeConstraintResidual(
        orbitals, work_orbitals, res, print_residual, norm_res);

    if (ct.isSpreadFunctionalEnergy()) addResidualSpreadPenalty(orbitals, res);

    comp_res_tm_.stop();

    return norm2Res;
}

//////////////////////////////////////////////////////////////////////////////
// compute res using psi and hpsi
template <class T>
void MGmol<T>::computeResidualUsingHPhi(T& psi,
    const T& hphi, T& res, const bool applyB)
{
    assert(psi.isCompatibleWith(hphi));
    assert(psi.isCompatibleWith(res));

    get_res_tm_.start();

    proj_matrices_->updateSubMatT();

    SquareLocalMatrices<MATDTYPE>& localT(proj_matrices_->getLocalT());

    pb::Lap<ORBDTYPE>* lapop = hamiltonian_->lapOper();
    const int ncolors        = psi.chromatic_number();
    if (ncolors > 0)
    {
        // compute B*psi and store in tmp
        ORBDTYPE* old_storage = 0;
        vector<ORBDTYPE> tmp;
        if (applyB)
        {
            const int ld = psi.getLda();
            tmp.resize(ld * ncolors);
            psi.setDataWithGhosts();
            psi.trade_boundaries();
            for (int i = 0; i < ncolors; i++)
            {
                lapop->rhs(psi.getFuncWithGhosts(i), &tmp[0] + ld * i);
            }

            // psi points to tmp temporarily
            old_storage = psi.getPsi(0);
            psi.set_storage(&tmp[0]);
        }

        // get B*phi*theta and store it in res in [Ry]
        // (even if Ritz functions mode)
        psi.multiplyByMatrix(localT, res);

        if (applyB)
        {
            psi.set_storage(old_storage);
        }

        // res = (B*phi*theta - H*phi) in [Ry]
        res.axpy(-1., hphi);
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

template <class T>
double MGmol<T>::computeConstraintResidual(T& orbitals,
    const T& hphi, T& res,
    const bool print_residual, const bool compute_norm_res)
{
    Control& ct(*(Control::instance()));

    computeResidualUsingHPhi(orbitals, hphi, res, ct.Mehrstellen());

    res.applyCorrMask();

    double normRes = -1.;
    if (print_residual || compute_norm_res)
    {
        if (ct.checkMaxResidual())
        {
            normRes = 0.5 * res.maxAbsValue(); // factor 0.5 for Rydberg to
                                               // Hartree conversion
            if (onpe0 && print_residual)
                os_ << setprecision(2) << scientific
                    << "max. |Residual_ij| = " << normRes << endl;
        }
        else
        {
            double norm2Res = res.dotProduct(res);
            normRes = 0.5 * sqrt(norm2Res); // factor 0.5 for Rydberg to Hartree
                                            // conversion
            if (onpe0 && print_residual)
                os_ << setprecision(2)
                    << scientific
                    //                   <<"|| Relative residual || =
                    //                   "<<normRes/ct.getNel()<<endl;
                    << "|| Residual || = " << normRes << endl;
        }
    }

    return normRes;
}

//////////////////////////////////////////////////////////////////////////////
// Get preconditioned residual in res_orbitals
template <class T>
double MGmol<T>::computePrecondResidual(T& phi,
    T& hphi, T& res, Ions& ions,
    KBPsiMatrixSparse* kbpsi, const bool print_residual, const bool norm_res)

{
    Control& ct = *(Control::instance());

    proj_matrices_->computeInvB();

    Potentials& pot          = hamiltonian_->potential();
    pb::Lap<ORBDTYPE>* lapop = hamiltonian_->lapOper();

    setGamma(*lapop, pot);

    // get H*psi stored in hphi
    // and psi^T H psi in Hij
    getHpsiAndTheta(ions, phi, hphi, kbpsi);

    double norm2Res
        = computeConstraintResidual(phi, hphi, res, print_residual, norm_res);

    if ((ct.getPrecondType() % 10) == 0 && ct.getMGlevels() >= 0)
    {
        // PRECONDITIONING
        // compute the preconditioned steepest descent direction
        // -> res
        orbitals_precond_->precond_mg(res);
    }

    // if( ct.isSpreadFunctionalActive() )addResidualSpreadPenalty(phi,res);

    return norm2Res;
}

//    Function to update potentials vh and vxc:
//
//    The new potentials are computed as a linear combination
//    of the old ones (input "vh" and "vxc") and the ones
//    corresponding to the input "rho".
template <class T>
void MGmol<T>::update_pot(const pb::GridFunc<POTDTYPE>& vh_init, const Ions& ions)
{
    electrostat_->setupInitialVh(vh_init);
    update_pot(ions);
}

template <class T>
void MGmol<T>::update_pot(const Ions& ions)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "Update potentials" << endl;
#endif

    Control& ct     = *(Control::instance());
    Potentials& pot = hamiltonian_->potential();

    // Update exchange-correlation potential
    xcongrid_->update();

    // Generate new hartree potential
    electrostat_->computeVh(ions, *rho_, pot);

    const bool flag_mixing = (fabs(ct.mix_pot - 1.) > 1.e-3);

    // evaluate potential correction
    if (flag_mixing)
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

template <class T>
void MGmol<T>::addResidualSpreadPenalty(T& phi, T& res)
{
    assert(spread_penalty_ != 0);

    spread_penalty_->addResidual(phi, res);
}

template class MGmol<LocGridOrbitals>;
template class MGmol<ExtendedGridOrbitals>;
