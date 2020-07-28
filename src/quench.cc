// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
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

#include "LocGridOrbitals.h"

#include "ABPG.h"
#include "AOMMprojector.h"
#include "Control.h"
#include "DFTsolver.h"
#include "DMStrategy.h"
#include "DavidsonSolver.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "KBPsiMatrixSparse.h"
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
Timer updateCenters_tm("MGmol<T>::updateCenters");

template <>
void MGmol<ExtendedGridOrbitals>::adaptLR(
    const SpreadsAndCenters<ExtendedGridOrbitals>* /*spreadf*/,
    const OrbitalsTransform* /*ot*/)
{
}

// depending on the value of ct.lr_updates_type, update
// localization regions:
// 0 -> center only
// 1 -> radius only
// 2 -> center and radius
template <class T>
void MGmol<T>::adaptLR(
    const SpreadsAndCenters<T>* spreadf, const OrbitalsTransform* ot)
{
    assert(lrs_);

    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp(" Adapt localization regions...", os_);

    if (spreadf == nullptr)
    {
        if (onpe0)
            os_ << "WARNING: Cannot adapt LRs: needs SpreadsAndCenters"
                << std::endl;
        return;
    }

    adaptLR_tm_.start();

    updateCenters_tm.start();
    // update centers
    if (ct.lr_updates_type != 1)
    {
        if (onpe0 && ct.verbose > 1)
            os_ << " Adapt localization centers..." << std::endl;
        lrs_->move(*spreadf);
    }
    updateCenters_tm.stop();
    // update radii
    if (ct.wannier_transform_type > 0 && ct.lr_updates_type > 0)
    {

        if (onpe0 && ct.verbose > 1)
            os_ << " Adapt localization radii..." << std::endl;
        const double vol_lrs = lrs_->computeVolume();
        if (onpe0 && ct.verbose > 1)
        {
            os_ << " current Volume LR = " << vol_lrs << std::endl;
        }

        double avg;
        // if ct.lr_volume_calc true, calculate volume base on spreads
        if (ct.lr_volume_calc == 1)
        {
            assert(ot != nullptr);
            const double vol_rotated = ot->volume();
            const double vol_spreads = spreadf->volume();
            if (onpe0 && ct.verbose > 1)
            {
                os_ << " vol. rotated localized orbitals = " << vol_rotated
                    << std::endl;
                os_ << " vol_spreads = " << vol_spreads << std::endl;
            }
            // target volume is vol_rotated*(ct.cut_radius**3)
            const double ratio
                = ct.cut_radius * cbrt(vol_rotated / vol_spreads);
            avg = lrs_->updateRadii(*spreadf, ratio);
        }
        else if (ct.lr_volume_calc == 2)
        {
            assert(ot != nullptr);
            if (onpe0) os_ << " Adapt LR using NOLMO spreads" << std::endl;
            const double ratio = ct.cut_radius;
            avg                = lrs_->updateRadii(ot, ratio);
        }
        else
        {
            if (onpe0) os_ << " Adapt with constant LR volume" << std::endl;
            avg = lrs_->updateRadiiConstVol(*spreadf);
        }

        if (onpe0 && ct.verbose > 1)
        {
            os_ << " Average radius = " << avg << std::endl;
        }
    }
    if (ct.verbose > 2) lrs_->printAllRegions(os_);

    if (ct.verbose > 1) lrs_->printMinMaxRadius(os_);

    adaptLR_tm_.stop();
}

template <class T>
void MGmol<T>::updateHmatrix(T& orbitals, Ions& ions)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "updateHmatrix()" << std::endl;
#endif

    // compute Hij
    getKBPsiAndHij(orbitals, ions);

    // save v[rho] used for H matrix to have a consistent energy
    // evaluation using H matrix
    energy_->saveVofRho();
}

template <class T>
void MGmol<T>::resetProjectedMatricesAndDM(T& orbitals, Ions& ions)
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
template <class T>
bool MGmol<T>::rotateStatesPairsCommonCenter(T& orbitals, T& work_orbitals)
{
    Control& ct        = *(Control::instance());
    bool wannier_pairs = false;

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    Vector3D ll(mygrid.ll(0), mygrid.ll(1), mygrid.ll(2));
    Vector3D origin(mygrid.origin(0), mygrid.origin(1), mygrid.origin(2));

    int st1, st2;
    std::map<int, Vector3D> save_centers;
    std::set<SymmetricPair> exclude_set;
    double drmin = lrs_->getStatesWithClosestCenters(&st1, &st2, exclude_set);
    if (onpe0)
        os_ << "drmin=" << drmin
            << ", ct.getMinDistanceCenters()=" << ct.getMinDistanceCenters()
            << std::endl;

    // if one distance smaller than tolerance, increase
    // tolerance to deal with more close centers right now
    // and avoid moving centers too often
    double tolerance = ct.getMinDistanceCenters();
    if (drmin < tolerance) tolerance = 1.2 * tolerance;

    while (drmin < tolerance)
    {
        if (onpe0 && ct.verbose > 1)
            os_ << "Min. distance between centers " << st1 << " and " << st2
                << " = " << drmin << std::endl;
        SpreadsAndCenters<T> spreadf2st(origin, ll);

        orbitals.orthonormalize2states(st1, st2);

        spreadf2st.computeSinCosDiag2states(orbitals, st1, st2);

        if (onpe0) spreadf2st.print(os_);

        // save centers
        save_centers.insert(
            std::pair<int, Vector3D>(st1, spreadf2st.computeCenter(0)));
        save_centers.insert(
            std::pair<int, Vector3D>(st2, spreadf2st.computeCenter(1)));

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
        for (std::map<int, Vector3D>::iterator it = save_centers.begin();
             it != save_centers.end(); it++)
            lrs_->setCenter(it->first, it->second);

        lrs_->setup();
    }

    return wannier_pairs;
}

// try to use some rotations to avoid degeneracies
template <class T>
bool MGmol<T>::rotateStatesPairsOverlap(
    T& orbitals, T& work_orbitals, const double threshold)
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
        os_ << "-------------------------" << std::endl;
        os_ << "rotateStatesPairsOverlap: compute smallest eigenvalue of S..."
            << std::endl;
    }
    double eigmin
        = proj_matrices_->getLinDependent2states(st1, st2, ct.verbose > 2);

    if (onpe0)
    {
        os_ << std::setprecision(3) << std::scientific << "eigmin=" << eigmin
            << ", threshold=" << threshold << std::endl;
    }

    while (eigmin < threshold)
    {
        if (onpe0 && ct.verbose > 1)
        {
            os_ << "-------------------------" << std::endl;
            os_ << std::setprecision(5) << std::scientific
                << "Min. eigenvalue for states " << st1 << ", " << st2 << " = "
                << eigmin << std::endl;
        }

        SpreadsAndCenters<T> spreadf2st(origin, ll);
        spreadf2st.computeSinCosDiag2states(orbitals, st1, st2);

        spreadf2st.print(os_);

        double d = spreadf2st.computeDistance(0, 1);
        if (d > ct.getThresholdDistancePairMLWF())
        {
            if (onpe0)
            {
                std::cout
                    << "Smallest eigenvalue not related to close LR centers..."
                    << std::endl;
                std::cout << "d=" << d << std::endl;
                std::cout << "stop MLWF transforms..." << std::endl;
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
            std::cout << "Gram Matrix for pair of states after transformation:"
                      << std::endl;
        if (projmatrices)
            projmatrices->printGramMatrix2states(st1, st2, std::cout);
        mmpi.barrier();
        if (onpe0 && ct.verbose > 1)
        {
            os_ << "-------------------------" << std::endl;
            os_ << "compute new smallest eigenvalue..." << std::endl;
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
                std::cout << "Same pair found twice in a row... stop MLWF "
                             "transforms..."
                          << std::endl;
            break;
        }
    };
    if (onpe0) std::cout << std::endl;

    return pairs;
}

template <class T>
void MGmol<T>::disentangleOrbitals(
    T& orbitals, T& work_orbitals, Ions& ions, int& max_steps)
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

template <>
void MGmol<LocGridOrbitals>::applyAOMMprojection(LocGridOrbitals& orbitals)
{
    aomm_ = new AOMMprojector(orbitals, lrs_);
    aomm_->projectOut(orbitals);
}

template <class T>
void MGmol<T>::applyAOMMprojection(T&)
{
}

template <>
int MGmol<LocGridOrbitals>::outerSolve(LocGridOrbitals& orbitals,
    LocGridOrbitals& work_orbitals, Ions& ions, const int max_steps,
    const int iprint, double& last_eks)
{
    int retval
        = 1; // 0 -> converged, -1 -> problem, -2 -> ( de>tol_energy_stop )
    Control& ct = *(Control::instance());
    // solve electronic structure problem
    // (inner iterations)
    switch (ct.OuterSolver())
    {
        case OuterSolverType::ABPG:
        case OuterSolverType::NLCG:
        {
            DFTsolver<LocGridOrbitals> solver(hamiltonian_, proj_matrices_,
                energy_, electrostat_, this, ions, rho_, dm_strategy_, os_);

            retval = solver.solve(
                orbitals, work_orbitals, ions, max_steps, iprint, last_eks);

            break;
        }

        case OuterSolverType::PolakRibiere:
        {
            PolakRibiereSolver<LocGridOrbitals> solver(hamiltonian_,
                proj_matrices_, energy_, electrostat_, this, ions, rho_,
                dm_strategy_, os_);

            retval = solver.solve(
                orbitals, work_orbitals, ions, max_steps, iprint, last_eks);

            break;
        }

        default:
            std::cerr << "Undefined iterative outer solver" << std::endl;
            return -1;
    }

    return retval;
}

template <class T>
int MGmol<T>::outerSolve(T& orbitals, T& work_orbitals, Ions& ions,
    const int max_steps, const int iprint, double& last_eks)
{
    int retval
        = 1; // 0 -> converged, -1 -> problem, -2 -> ( de>tol_energy_stop )
    Control& ct = *(Control::instance());
    // solve electronic structure problem
    // (inner iterations)
    switch (ct.OuterSolver())
    {
        case OuterSolverType::ABPG:
        case OuterSolverType::NLCG:
        {
            DFTsolver<T> solver(hamiltonian_, proj_matrices_, energy_,
                electrostat_, this, ions, rho_, dm_strategy_, os_);

            retval = solver.solve(
                orbitals, work_orbitals, ions, max_steps, iprint, last_eks);

            break;
        }

        case OuterSolverType::PolakRibiere:
        {
            PolakRibiereSolver<T> solver(hamiltonian_, proj_matrices_, energy_,
                electrostat_, this, ions, rho_, dm_strategy_, os_);

            retval = solver.solve(
                orbitals, work_orbitals, ions, max_steps, iprint, last_eks);

            break;
        }

        case OuterSolverType::Davidson:
        {
            const std::vector<std::vector<int>>& gids(
                orbitals.getOverlappingGids());

            DavidsonSolver<T, dist_matrix::DistMatrix<DISTMATDTYPE>> solver(
                comm_, os_, *ions_, hamiltonian_, rho_, energy_, electrostat_,
                this, ct.numst, ct.occ_width, ct.getNel(), gids);

            retval = solver.solve(orbitals, work_orbitals);
            break;
        }

        default:
            std::cerr << "Undefined iterative outer solver" << std::endl;
            return -1;
    }

    return retval;
}

template <class T>
int MGmol<T>::quench(T* orbitals, Ions& ions, const int max_inner_steps,
    const int iprint, double& last_eks)
{
    assert(max_inner_steps > -1);

    quench_tm.start();
    int retval
        = 1; // 0 -> converged, -1 -> problem, -2 -> ( de>tol_energy_stop )
    int max_steps = max_inner_steps;

    Control& ct = *(Control::instance());
    if (ct.restart_info > 2)
    {
        DFTsolver<T>::setItCountLarge();
    }

    if (onpe0)
    {
        os_ << "###########################" << std::endl;
        os_ << "QUENCH ELECTRONS" << std::endl;
        os_ << "###########################" << std::endl;
    }

    // get actual indexes of stored functions
    const std::vector<std::vector<int>>& gids(orbitals->getOverlappingGids());

    g_kbpsi_->setup(*ions_);
    electrostat_->setup(ct.vh_its);
    rho_->setup(ct.getOrbitalsType(), gids);

    T work_orbitals("Work", *orbitals);

    orbitals->setDataWithGhosts();
    orbitals->trade_boundaries();

    disentangleOrbitals(*orbitals, work_orbitals, ions, max_steps);

    // setup "kernel" functions for AOMM algorithm
    if (ct.use_kernel_functions)
    {
        applyAOMMprojection(*orbitals);
    }

    orbitals_precond_ = new OrbitalsPreconditioning<T>();
    orbitals_precond_->setup(
        *orbitals, ct.getMGlevels(), ct.lap_type, currentMasks_, lrs_);

    // solve electronic structure problem
    // (inner iterations)
    retval = outerSolve(
        *orbitals, work_orbitals, ions, max_steps, iprint, last_eks);
    if (retval == -1) return -1;

    if (ct.use_kernel_functions)
    {
        delete aomm_;
        aomm_ = nullptr;
    }
    delete orbitals_precond_;
    orbitals_precond_ = nullptr;

    // Get the n.l. energy
    // TODO: Fix bug where energy vs. time output is incorrect if get_evnl is
    // not called.
    // if( ct.verbose>1 && ct.short_sighted )
    if (ct.short_sighted)
    {
        const double evnl = get_evnl(ions, *orbitals);
        if (onpe0)
        {
            os_ << std::setprecision(8) << std::fixed
                << " Enl            [Ha] = " << evnl << std::endl;
        }
    }

    const double ts = 0.5 * proj_matrices_->computeEntropy(); // in [Ha]
    if (ct.verbose > 1)
    {
        if (onpe0)
        {
            os_ << std::setprecision(8) << std::fixed
                << " TS             [Ha] = " << ts << std::endl;
        }
    }
    last_eks = energy_->evaluateTotal(ts, proj_matrices_, *orbitals, 2, os_);

    if (ct.computeCondGramMD())
    {
        double condS = proj_matrices_->computeCond();
        if (onpe0)
            os_ << std::setprecision(2) << std::scientific
                << "Condition Number of S: " << condS << std::endl;
    }

    if (ct.isLocMode() || ct.isSpreadFunctionalActive())
    {
        // build matrices necessary to compute spreads and centers
        spreadf_->computePositionMatrix(*orbitals, work_orbitals);

        if (ct.verbose > 0 || !ct.AtomsMove())
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

template class MGmol<LocGridOrbitals>;
template class MGmol<ExtendedGridOrbitals>;
