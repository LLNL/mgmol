// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "MVPSolver.h"
#include "Control.h"
#include "DistMatrixWithSparseComponent.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "Ions.h"
#include "KBPsiMatrixSparse.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "Potentials.h"
#include "ProjectedMatrices.h"
#include "Rho.h"
#include "tools.h"

#include <iomanip>

Timer MVPSolver::solve_tm_("MVPSolver::solve");
Timer MVPSolver::target_tm_("MVPSolver::target");

static int sparse_distmatrix_nb_partitions = 128;

double evalEntropyMVP(ProjectedMatricesInterface* projmatrices,
    const bool print_flag, ostream& os)
{
    const double ts = 0.5 * projmatrices->computeEntropy(); // in [Ha]
    if (onpe0 && print_flag)
        os << setprecision(8) << "compute Entropy: TS=" << ts << "[Ha]" << endl;

    return ts;
}

MVPSolver::MVPSolver(MPI_Comm comm, ostream& os, Ions& ions, Rho* rho,
    Energy* energy, Electrostatic* electrostat, MGmol* mgmol_strategy,
    const int numst, const double kbT, const int nel,
    const vector<vector<int>>& global_indexes, const short n_inner_steps,
    const bool use_old_dm)
    : comm_(comm),
      os_(os),
      ions_(ions),
      n_inner_steps_(n_inner_steps),
      use_old_dm_(use_old_dm)
{
    history_length_ = 2;
    eks_history_.resize(history_length_, 100000.);

    rho_            = rho;
    energy_         = energy;
    electrostat_    = electrostat;
    mgmol_strategy_ = mgmol_strategy;

    numst_ = numst;
    work_
        = new dist_matrix::DistMatrix<DISTMATDTYPE>("workMVP", numst_, numst_);

    proj_mat_work_ = new ProjectedMatrices(numst_, false);
    proj_mat_work_->setup(kbT, nel, global_indexes);
}

MVPSolver::~MVPSolver()
{
    delete work_;
    delete proj_mat_work_;
}

double MVPSolver::evaluateDerivative(
    dist_matrix::DistMatrix<DISTMATDTYPE>& dmInit,
    dist_matrix::DistMatrix<DISTMATDTYPE>& delta_dm, const double ts0)
{
    work_->symm('l', 'l', 1., proj_mat_work_->getMatHB(), delta_dm, 0.);

    // factor 0.5 to get value in Hartree (HB is in Rydberg)
    double de0  = 0.5 * work_->trace();
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
        os_ << "Derivative of U in 0 = " << de0 << endl;

    //
    // evaluate numerical derivative of entropy in beta=0
    //
    // if( onpe0 ) os_<<"evaluate numerical derivative of entropy in
    // beta=0"<<endl;
    const double dbeta = 0.001;
    *work_             = dmInit;
    work_->axpy(dbeta, delta_dm);
    // proj_mat_work_->setDM(*work_,orbitals.getIterativeIndex());
    proj_mat_work_->setDM(*work_, -1);
    proj_mat_work_->computeOccupationsFromDM();

    const double tsd0e = evalEntropyMVP(proj_mat_work_, false, os_);
    const double dts0e = (tsd0e - ts0) / dbeta;

    // if( onpe0 )os_<<"ts0="<<scientific<<ts0<<endl;
    // if( onpe0 )os_<<"tsd0e with 1.*beta="<<scientific<<tsd0e<<endl;

    if (onpe0 && ct.verbose > 2)
    {
        os_ << setprecision(8);
        os_ << "Numerical derivative of -TS in 0 = " << scientific << -dts0e
            << endl;
        // os_<<"interval for finite differences = "<<dbeta<<endl;
    }
    de0 -= dts0e;

    return de0;
}

void MVPSolver::buildTarget_MVP(dist_matrix::DistMatrix<DISTMATDTYPE>& h11,
    dist_matrix::DistMatrix<DISTMATDTYPE>& s11,
    dist_matrix::DistMatrix<DISTMATDTYPE>& target)
// build target
{
    target_tm_.start();

    Control& ct = *(Control::instance());

    proj_mat_work_->assignH(h11);

    int orbitals_index = 0;
    proj_mat_work_->setGramMatrix(s11, orbitals_index);

    //
    // compute target density matrix, with occupations>0 in numst only
    //
    // if( onpe0 )os_<<"Compute XN..."<<endl;
    proj_mat_work_->setHB2H();

    // if( onpe0 )os_<<"Compute F with Nel="<<ct.getNel()<<endl;
    proj_mat_work_->updateDM(orbitals_index);

    target = proj_mat_work_->dm();

    if (ct.verbose > 2)
    {
        double pnel = proj_mat_work_->getNel();
        if (onpe0) os_ << "Number of electrons = " << pnel << endl;
    }

    target_tm_.stop();
}

// update density matrix in N x N space
int MVPSolver::solve(LocGridOrbitals& orbitals)
{
    Control& ct = *(Control::instance());

    assert(numst_ = (int)orbitals.numst());
    assert(n_inner_steps_ > 0);

    solve_tm_.start();

    if (onpe0 && ct.verbose > 1)
    {
        os_ << "---------------------------------------------------------------"
               "-"
            << endl;
        os_ << "Update DM functions using MVP Solver..." << endl;
        os_ << "---------------------------------------------------------------"
               "-"
            << endl;
    }

    // step for numerical derivative

    // double nel=orbitals.projMatrices()->getNel();
    // if( onpe0 )
    //    os_<<"NEW algorithm: Number of electrons at start = "<<nel<<endl;

    KBPsiMatrixSparse kbpsi(0);
    kbpsi.setup(ions_, orbitals);

    dist_matrix::DistMatrix<DISTMATDTYPE> s11("s11", numst_, numst_);
    {
        ProjectedMatrices* projmatrices
            = dynamic_cast<ProjectedMatrices*>(orbitals.getProjMatrices());
        s11 = projmatrices->getGramMatrix();
    }

    {
        dist_matrix::DistMatrix<DISTMATDTYPE> dmInit("dm", numst_, numst_);

        // save computed vh for a fair energy "comparison" with vh computed
        // in close neigborhood
        const pb::GridFunc<POTDTYPE> vh_init(electrostat_->getVh());

        orbitals.setDataWithGhosts();

        // Update density
        if (use_old_dm_) rho_->update(orbitals);

        // compute linear component of H
        dist_matrix::DistMatrixWithSparseComponent<DISTMATDTYPE> h11("h11",
            numst_, numst_, comm_, mgmol_strategy_->getRemoteTasksDistMatrix(),
            sparse_distmatrix_nb_partitions);

        kbpsi.computeAll(ions_, orbitals);

        kbpsi.computeHvnlMatrix(&kbpsi, ions_, h11);

        for (int inner_it = 0; inner_it < n_inner_steps_; inner_it++)
        {
            if (onpe0 && ct.verbose > 1)
            {
                os_ << "---------------------------" << endl;
                os_ << "Inner iteration " << inner_it << endl;
                os_ << "---------------------------" << endl;
            }

            // Update potential
            if (use_old_dm_ || inner_it > 0)
            {
                mgmol_strategy_->update_pot(vh_init, ions_);

                energy_->saveVofRho();
            }

            ProjectedMatricesInterface* current_proj_mat
                = (inner_it == 0) ? orbitals.getProjMatrices() : proj_mat_work_;

            double ts0;
            double e0        = 0.;
            const int printE = (ct.verbose > 1) ? 1 : 0;

            // compute h11 for the current potential by adding local part to
            // nonlocal components (old local part reset to zero)
            mgmol_strategy_->addHlocal2matrix(orbitals, orbitals, h11);

            if (inner_it == 0)
            {
                // orbitals are new, so a few things need to recomputed
                ProjectedMatrices* projmatrices
                    = dynamic_cast<ProjectedMatrices*>(
                        orbitals.getProjMatrices());

                projmatrices->assignH(h11);
                projmatrices->setHB2H();

                if (use_old_dm_) dmInit = projmatrices->dm();

                ts0 = evalEntropyMVP(projmatrices, true, os_);
                e0  = energy_->evaluateTotal(
                    ts0, projmatrices, orbitals, printE, os_);
            }
            else
            {

                dmInit = proj_mat_work_->dm();

                proj_mat_work_->assignH(h11);
                proj_mat_work_->setHB2H();

                ts0 = evalEntropyMVP(proj_mat_work_, true, os_);
                e0  = energy_->evaluateTotal(
                    ts0, proj_mat_work_, orbitals, printE, os_);
            }

            // N x N target...
            // proj_mat_work_->setHiterativeIndex(orbitals.getIterativeIndex(),
            // pot.getIterativeIndex());

            dist_matrix::DistMatrix<DISTMATDTYPE> target(
                "target", numst_, numst_);
            dist_matrix::DistMatrix<DISTMATDTYPE> delta_dm(
                "delta_dm", numst_, numst_);

            buildTarget_MVP(h11, s11, target);

            if (use_old_dm_ || inner_it > 0)
            {
                delta_dm = target;
                delta_dm -= dmInit;

                double de0 = evaluateDerivative(dmInit, delta_dm, ts0);

                //
                // evaluate free energy at beta=1
                //
                if (onpe0 && ct.verbose > 2) cout << "Target energy..." << endl;
                proj_mat_work_->setDM(target, orbitals.getIterativeIndex());
                proj_mat_work_->computeOccupationsFromDM();
                if (ct.verbose > 2) current_proj_mat->printOccupations(os_);
                double nel = proj_mat_work_->getNel();
                if (onpe0 && ct.verbose > 1)
                    os_ << "Number of electrons at beta=1 : " << nel << endl;

                // if( onpe0 )os_<<"Rho..."<<endl;
                rho_->computeRho(orbitals, target);
                rho_->rescaleTotalCharge();

                mgmol_strategy_->update_pot(vh_init, ions_);

                energy_->saveVofRho();

                // update h11
                {
                    mgmol_strategy_->addHlocal2matrix(orbitals, orbitals, h11);
                }

                proj_mat_work_->assignH(h11);
                proj_mat_work_->setHB2H();

                const double ts1
                    = evalEntropyMVP(proj_mat_work_, (ct.verbose > 2), os_);
                const double e1 = energy_->evaluateTotal(
                    ts1, proj_mat_work_, orbitals, ct.verbose - 1, os_);

                // line minimization
                double beta
                    = minQuadPolynomial(e0, e1, de0, (ct.verbose > 2), os_);

                if (onpe0 && ct.verbose > 0)
                {
                    os_ << setprecision(12);
                    os_ << fixed << "Inner iteration " << inner_it
                        << ", E0=" << e0 << ", E1=" << e1;
                    os_ << scientific << ", E0'=" << de0 << " -> beta=" << beta;
                    os_ << endl;
                }

                // update DM
                *work_ = dmInit;
                work_->axpy(beta, delta_dm);
            }
            else // no old DM to use
            {
                *work_ = target;
            }
            proj_mat_work_->setDM(*work_, orbitals.getIterativeIndex());

            if (inner_it < n_inner_steps_ - 1)
            {
                if (ct.verbose > 2)
                {
                    double pnel = proj_mat_work_->getNel();
                    if (onpe0)
                        os_ << "Number of electrons for interpolated DM = "
                            << pnel << endl;
                }

                // if( onpe0 )os_<<"Rho..."<<endl;
                rho_->computeRho(orbitals, *work_);
                rho_->rescaleTotalCharge();
            }

        } // inner iterations

        ProjectedMatrices* projmatrices
            = dynamic_cast<ProjectedMatrices*>(orbitals.getProjMatrices());
        projmatrices->setDM(*work_, orbitals.getIterativeIndex());
    }

    // Generate new density
    rho_->update(orbitals);

    if (onpe0 && ct.verbose > 1)
    {
        os_ << "---------------------------------------------------------------"
               "-"
            << endl;
        os_ << "End MVP Solver..." << endl;
        os_ << "---------------------------------------------------------------"
               "-"
            << endl;
    }
    solve_tm_.stop();

    return 0;
}

void MVPSolver::printTimers(ostream& os)
{
    if (onpe0)
    {
        os << setprecision(2) << fixed << endl;
        solve_tm_.print(os);
        target_tm_.print(os);
    }
}
