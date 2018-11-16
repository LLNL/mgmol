// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <time.h>
#include <vector>
using namespace std;

#include "LocGridOrbitals.h"

#include "ABPG.h"
#include "AOMMprojector.h"
#include "Control.h"
#include "DFTsolver.h"
#include "DMStrategy.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "GrassmanCG.h"
#include "GrassmanCGSparse.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "KBPsiMatrix.h"
#include "LocalizationRegions.h"
#include "MGmol.h"
#include "MPIdata.h"
#include "MasksSet.h"
#include "Mesh.h"
#include "OrbitalsPreconditioning.h"
#include "OrbitalsTransform.h"
#include "PolakRibiereSolver.h"
#include "Potentials.h"
#include "ProjectedMatricesInterface.h"
#include "ProjectedMatricesSparse.h"
#include "Rho.h"
#include "SpreadsAndCenters.h"
#include "SubspaceProjector.h"
#include "SymmetricPair.h"
#include "tools.h"

#define TEST_ENERGY 0

Timer quench_tm("quench");
Timer quench_evnl_tm("quench_evnl");
Timer MGmol::adaptLR_tm_("MGmol::adaptLR");
Timer updateCenters_tm("MGmol::updateCenters");

// depending on the value of ct.lr_updates_type, update
// localization regions:
// 0 -> center only
// 1 -> radius only
// 2 -> center and radius
void MGmol::adaptLR(
    const SpreadsAndCenters* spreadf, const OrbitalsTransform* ot)
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp(" Adapt localization regions...", os_);

    if (spreadf == NULL)
    {
        if (onpe0)
            os_ << "WARNING: Cannot adapt LRs: needs SpreadsAndCenters" << endl;
        return;
    }

    adaptLR_tm_.start();

    updateCenters_tm.start();
    // update centers
    if (ct.lr_updates_type != 1)
    {
        if (onpe0 && ct.verbose > 1)
            os_ << " Adapt localization centers..." << endl;
        lrs_->move(*spreadf);
    }
    updateCenters_tm.stop();
    // update radii
    if (ct.wannier_transform_type > 0 && ct.lr_updates_type > 0)
    {

        if (onpe0 && ct.verbose > 1)
            os_ << " Adapt localization radii..." << endl;
        const double vol_lrs = lrs_->computeVolume();
        if (onpe0 && ct.verbose > 1)
        {
            os_ << " current Volume LR = " << vol_lrs << endl;
        }

        double avg;
        // if ct.lr_volume_calc true, calculate volume base on spreads
        if (ct.lr_volume_calc == 1)
        {
            assert(ot != NULL);
            const double vol_rotated = ot->volume();
            const double vol_spreads = spreadf->volume();
            if (onpe0 && ct.verbose > 1)
            {
                os_ << " vol. rotated localized orbitals = " << vol_rotated
                    << endl;
                os_ << " vol_spreads = " << vol_spreads << endl;
            }
            // target volume is vol_rotated*(ct.cut_radius**3)
            const double ratio
                = ct.cut_radius * cbrt(vol_rotated / vol_spreads);
            avg = lrs_->updateRadii(*spreadf, ratio);
        }
        else if (ct.lr_volume_calc == 2)
        {
            assert(ot != NULL);
            if (onpe0) os_ << " Adapt LR using NOLMO spreads" << endl;
            const double ratio = ct.cut_radius;
            avg                = lrs_->updateRadii(ot, ratio);
        }
        else
        {
            if (onpe0) os_ << " Adapt with constant LR volume" << endl;
            avg = lrs_->updateRadiiConstVol(*spreadf);
        }

        if (onpe0 && ct.verbose > 1)
        {
            os_ << " Average radius = " << avg << endl;
        }
    }
    if (ct.verbose > 2) lrs_->printAllRegions(os_);

    if (ct.verbose > 1) lrs_->printMinMaxRadius(os_);

    adaptLR_tm_.stop();
}

void MGmol::updateHmatrix(LocGridOrbitals& orbitals, Ions& ions)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "updateHmatrix()" << endl;
#endif

    // compute Hij
    getKBPsiAndHij(orbitals, ions);

    // save v[rho] used for H matrix to have a consistent energy
    // evaluation using H matrix
    energy_->saveVofRho();
}

void MGmol::resetProjectedMatricesAndDM(LocGridOrbitals& orbitals, Ions& ions)
{
    orbitals.computeGramAndInvS();

    // compute new H
    updateHmatrix(orbitals, ions);

    orbitals.computeBAndInvB(*hamiltonian_->lapOper());

    proj_matrices_->updateThetaAndHB();

    // reset DM using chosen strategy
    dm_strategy_->initialize();
}

// try to use some rotations to avoid degeneracies
bool MGmol::rotateStatesPairsCommonCenter(
    LocGridOrbitals& orbitals, LocGridOrbitals& work_orbitals)
{
    Control& ct        = *(Control::instance());
    bool wannier_pairs = false;

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    Vector3D ll(mygrid.ll(0), mygrid.ll(1), mygrid.ll(2));
    Vector3D origin(mygrid.origin(0), mygrid.origin(1), mygrid.origin(2));

    int st1, st2;
    map<int, Vector3D> save_centers;
    set<SymmetricPair> exclude_set;
    double drmin = lrs_->getStatesWithClosestCenters(&st1, &st2, exclude_set);
    if (onpe0)
        os_ << "drmin=" << drmin
            << ", ct.getMinDistanceCenters()=" << ct.getMinDistanceCenters()
            << endl;

    // if one distance smaller than tolerance, increase
    // tolerance to deal with more close centers right now
    // and avoid moving centers too often
    double tolerance = ct.getMinDistanceCenters();
    if (drmin < tolerance) tolerance = 1.2 * tolerance;

    while (drmin < tolerance)
    {
        if (onpe0 && ct.verbose > 1)
            os_ << "Min. distance between centers " << st1 << " and " << st2
                << " = " << drmin << endl;
        SpreadsAndCenters spreadf2st(origin, ll);

        orbitals.orthonormalize2states(st1, st2);

        // double n1=orbitals.normState(st1);
        // double n2=orbitals.normState(st2);
        // if( onpe0 && ct.verbose>1 )
        //    os_<<"n1="<<n1<<", n2="<<n2<<endl;

        spreadf2st.computeSinCosDiag2states(orbitals, st1, st2);

        if (onpe0) spreadf2st.print(os_);

        // save centers
        save_centers.insert(
            pair<int, Vector3D>(st1, spreadf2st.computeCenter(0)));
        save_centers.insert(
            pair<int, Vector3D>(st2, spreadf2st.computeCenter(1)));

        getMLWF2states(st1, st2, orbitals, work_orbitals);
        spreadf2st.computeSinCosDiag2states(orbitals, st1, st2);

        if (onpe0) spreadf2st.print(os_);

        SymmetricPair pair(st1, st2);
        exclude_set.insert(pair);

        // reset centers in lrs_
        // so that they are not found again in next step
        lrs_->setCenter(st1, spreadf2st.computeCenter(0));
        lrs_->setCenter(st2, spreadf2st.computeCenter(1));

        drmin = lrs_->getStatesWithClosestCenters(&st1, &st2, exclude_set);

        wannier_pairs = true;
    }

    if (wannier_pairs)
    {
        // now put back old centers into lrs_
        // to avoid inconsistencies with LR allocations
        for (map<int, Vector3D>::iterator it = save_centers.begin();
             it != save_centers.end(); it++)
            lrs_->setCenter(it->first, it->second);

        lrs_->setup();
    }

    return wannier_pairs;
}

// try to use some rotations to avoid degeneracies
bool MGmol::rotateStatesPairsOverlap(LocGridOrbitals& orbitals,
    LocGridOrbitals& work_orbitals, const double threshold)
{
    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    bool pairs      = false;

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    Vector3D ll(mygrid.ll(0), mygrid.ll(1), mygrid.ll(2));
    Vector3D origin(mygrid.origin(0), mygrid.origin(1), mygrid.origin(2));

    int st1, st2;

    if (onpe0 && ct.verbose > 2)
    {
        os_ << "-------------------------" << endl;
        os_ << "rotateStatesPairsOverlap: compute smallest eigenvalue of S..."
            << endl;
    }
    double eigmin
        = proj_matrices_->getLinDependent2states(st1, st2, ct.verbose > 2);

    if (onpe0)
    {
        os_ << setprecision(3) << scientific << "eigmin=" << eigmin
            << ", threshold=" << threshold << endl;
    }

    while (eigmin < threshold)
    {
        if (onpe0 && ct.verbose > 1)
        {
            os_ << "-------------------------" << endl;
            os_ << setprecision(5) << scientific
                << "Min. eigenvalue for states " << st1 << ", " << st2 << " = "
                << eigmin << endl;
        }

        SpreadsAndCenters spreadf2st(origin, ll);
        spreadf2st.computeSinCosDiag2states(orbitals, st1, st2);

        spreadf2st.print(os_);

        double d = spreadf2st.computeDistance(0, 1);
        if (d > ct.getThresholdDistancePairMLWF())
        {
            if (onpe0)
            {
                cout << "Smallest eigenvalue not related to close LR centers..."
                     << endl;
                cout << "d=" << d << endl;
                cout << "stop MLWF transforms..." << endl;
            }
            break;
        }

        getMLWF2states(st1, st2, orbitals, work_orbitals);
        spreadf2st.computeSinCosDiag2states(orbitals, st1, st2);

        spreadf2st.print(os_);
        orbitals.computeGramAndInvS();

        ProjectedMatricesSparse* projmatrices
            = dynamic_cast<ProjectedMatricesSparse*>(proj_matrices_);
        if (onpe0)
            cout << "Gram Matrix for pair of states after transformation:"
                 << endl;
        if (projmatrices) projmatrices->printGramMatrix2states(st1, st2, cout);
        mmpi.barrier();
        if (onpe0 && ct.verbose > 1)
        {
            os_ << "-------------------------" << endl;
            os_ << "compute new smallest eigenvalue..." << endl;
        }

        const int oldst1 = st1;
        const int oldst2 = st2;
        // compute new smallest eigenvalue
        eigmin
            = proj_matrices_->getLinDependent2states(st1, st2, ct.verbose > 2);
        // remove 5% to deal with other small eigenvalues while we are at it
        // (perturbing solution...)
        eigmin *= 0.95;

        pairs = true;
        if (oldst1 == st1 && oldst2 == st2)
        {
            if (onpe0)
                cout << "Same pair found twice in a row... stop MLWF "
                        "transforms..."
                     << endl;
            break;
        }
    };
    if (onpe0) cout << endl;

    return pairs;
}

void MGmol::disentangleOrbitals(LocGridOrbitals& orbitals,
    LocGridOrbitals& work_orbitals, Ions& ions, int& max_steps)
{
    Control& ct = *(Control::instance());

    bool wannier_pairs = false;

    // rotate pairs if smallest eigenvalue of overlap matrix below threshold
    if (ct.getThresholdEigenvalueGram() > 0.)
    {
        wannier_pairs = rotateStatesPairsOverlap(
            orbitals, work_orbitals, ct.getThresholdEigenvalueGram());
    }
    if (ct.getMinDistanceCenters() > 0. && ct.isLocMode())
    {
        wannier_pairs = rotateStatesPairsCommonCenter(orbitals, work_orbitals);
    }

    if (wannier_pairs)
    {
        // clear old centers to avoid extrapolation with non-matching centers
        lrs_->clearOldCenters();

        resetProjectedMatricesAndDM(orbitals, ions);
        // increase number of inner steps to deal with local reorthogonalization
        max_steps *= 2;
        ct.num_MD_steps--;
    }
}

int MGmol::quench(LocGridOrbitals* orbitals, Ions& ions,
    const int max_inner_steps, const int iprint, double& last_eks)
{
    assert(max_inner_steps > -1);

    quench_tm.start();
    int retval
        = 1; // 0 -> converged, -1 -> problem, -2 -> ( de>tol_energy_stop )
    int max_steps = max_inner_steps;

    Control& ct = *(Control::instance());
    if (ct.restart_info > 2)
    {
        DFTsolver::setItCountLarge();
    }

    if (onpe0)
    {
        os_ << "###########################" << endl;
        os_ << "QUENCH ELECTRONS" << endl;
        os_ << "###########################" << endl;
    }

    // get actual indexes of stored functions
    const vector<vector<int>>& gids(orbitals->getOverlappingGids());

    g_kbpsi_->setup(*ions_, *orbitals);
    electrostat_->setup(ct.vh_its);
    rho_->setup(ct.orbital_type, gids);

    LocGridOrbitals work_orbitals(*orbitals);

    orbitals->setDataWithGhosts();
    orbitals->trade_boundaries();

    disentangleOrbitals(*orbitals, work_orbitals, ions, max_steps);

    // setup "kernel" functions for AOMM algorithm
    if (ct.use_kernel_functions)
    {
        aomm_ = new AOMMprojector(*orbitals, *lrs_);
        aomm_->projectOut(*orbitals);
    }

    orbitals_precond_ = new OrbitalsPreconditioning();
    orbitals_precond_->setup(
        *orbitals, ct.getMGlevels(), ct.lap_type, currentMasks_, lrs_);

    // solve electronic structure problem
    // (inner iterations)
    switch (ct.it_algo_type)
    {
        case 0:
        {
            DFTsolver solver(hamiltonian_, proj_matrices_, energy_,
                electrostat_, this, ions, rho_, dm_strategy_, os_);

            retval = solver.solve(
                *orbitals, work_orbitals, ions, max_steps, iprint, last_eks);

            break;
        }

        case 1:
        {
            DFTsolver solver(hamiltonian_, proj_matrices_, energy_,
                electrostat_, this, ions, rho_, dm_strategy_, os_);

            retval = solver.solve(
                *orbitals, work_orbitals, ions, max_steps, iprint, last_eks);

            break;
        }

        case 3:
        {
            PolakRibiereSolver solver(hamiltonian_, proj_matrices_, energy_,
                electrostat_, this, ions, rho_, dm_strategy_, os_);

            retval = solver.solve(
                *orbitals, work_orbitals, ions, max_steps, iprint, last_eks);

            break;
        }

        default:
            cerr << "Undefined iterative electronic structure solver" << endl;
            return -1;
    }

    if (ct.use_kernel_functions)
    {
        delete aomm_;
        aomm_ = 0;
    }
    delete orbitals_precond_;
    orbitals_precond_ = 0;

    // Get the n.l. energy
    // TODO: Fix bug where energy vs. time output is incorrect if get_evnl is
    // not called.
    // if( ct.verbose>1 && ct.short_sighted )
    if (ct.short_sighted)
    {
        const double evnl = get_evnl(ions, *orbitals);
        if (onpe0)
        {
            os_ << setprecision(8) << fixed << " Enl            [Ha] = " << evnl
                << endl;
        }
    }

    const double ts = 0.5 * proj_matrices_->computeEntropy(); // in [Ha]
    if (ct.verbose > 1)
    {
        if (onpe0)
        {
            os_ << setprecision(8) << fixed << " TS             [Ha] = " << ts
                << endl;
        }
    }
    last_eks = energy_->evaluateTotal(ts, proj_matrices_, *orbitals, 2, os_);

    if (ct.computeCondGramMD())
    {
        double condS = proj_matrices_->computeCond();
        if (onpe0)
            os_ << setprecision(2) << scientific
                << "Condition Number of S: " << condS << endl;
    }

    if (ct.isLocMode() || ct.isSpreadFunctionalActive())
    {
        // build matrices necessary to compute spreads and centers
        spreadf_->computePositionMatrix(*orbitals, work_orbitals);

        if (ct.verbose > 0 || ct.atoms_dyn == 0)
        {
            spreadf_->print(os_);
            spreadf_->printStats(os_);
        }
    }

    if (!ct.isLocMode())
    {
        if (ct.wannier_transform_type >= 1)
        {
            wftransform(orbitals, &work_orbitals, ions);
        }
    }

    // delete hnl_orbitals_;hnl_orbitals_=0;

    quench_tm.stop();
    return retval;
}
