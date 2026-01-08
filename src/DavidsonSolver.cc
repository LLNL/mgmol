// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// the Lawrence Livermore National Laboratory.
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DavidsonSolver.h"
#include "Control.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "ExtendedGridOrbitals.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "KBPsiMatrixSparse.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "Potentials.h"
#include "ProjectedMatrices2N.h"
#include "ReplicatedMatrix.h"
#include "Rho.h"
#include "tools.h"

#include <iomanip>

template <class OrbitalsType, class MatrixType>
Timer DavidsonSolver<OrbitalsType, MatrixType>::solve_tm_(
    "DavidsonSolver::solve");
template <class OrbitalsType, class MatrixType>
Timer DavidsonSolver<OrbitalsType, MatrixType>::target_tm_(
    "DavidsonSolver::target");

double evalEntropy(ProjectedMatricesInterface* projmatrices,
    const bool print_flag, std::ostream& os)
{
    const double ts = 0.5 * projmatrices->computeEntropy(); // in [Ha]
    if (onpe0 && print_flag)
        os << std::setprecision(8) << "compute Entropy: TS=" << ts << "[Ha]"
           << std::endl;

    return ts;
}

template <class OrbitalsType, class MatrixType>
DavidsonSolver<OrbitalsType, MatrixType>::DavidsonSolver(std::ostream& os,
    Ions& ions, Hamiltonian<OrbitalsType>* hamiltonian, Rho<OrbitalsType>* rho,
    Energy<OrbitalsType>* energy, Electrostatic* electrostat,
    MGmol<OrbitalsType>* mgmol_strategy,
    const std::vector<std::vector<int>>& global_indexes, const double mixing,
    const bool with_spin)
    : os_(os),
      ions_(ions),
      hamiltonian_(hamiltonian),
      rho_(rho),
      energy_(energy),
      electrostat_(electrostat),
      mgmol_strategy_(mgmol_strategy),
      mixing_(mixing),
      diel_control_(hamiltonian->potential(), electrostat, os, onpe0)
{
    Control& ct = *(Control::instance());

    history_length_ = 2;
    eks_history_.resize(history_length_, 100000.);

    numst_ = ct.numst;
    work2N_.reset(new MatrixType("work2N", 2 * numst_, 2 * numst_));

    proj_mat2N_.reset(new ProjectedMatrices2N<MatrixType>(
        2 * ct.numst, with_spin, ct.occ_width));
    proj_mat2N_->setup(global_indexes);
}

template <class OrbitalsType, class MatrixType>
DavidsonSolver<OrbitalsType, MatrixType>::~DavidsonSolver()
{
}

template <class OrbitalsType, class MatrixType>
int DavidsonSolver<OrbitalsType, MatrixType>::checkConvergence(
    const double e0, const int it, const double tol)
{
    // save energy recent history
    for (int i = history_length_ - 1; i > 0; i--)
        eks_history_[i] = eks_history_[i - 1];

    eks_history_[0] = e0;
    de_old_         = de_;
    de_             = fabs(eks_history_[0] - eks_history_[1]);
    double de2      = std::max(de_, de_old_);
    if (onpe0)
    {
        os_ << std::setprecision(12) << std::fixed << "%%    " << it
            << " SC ENERGY = " << eks_history_[0];
        os_ << std::setprecision(2) << std::scientific
            << ", delta Eks = " << de_ << std::endl;
    }

    // test if this iteration seems converged
    int retval = 1;
    if (it > 2 && de2 < tol)
    {
        if (onpe0)
            os_ << std::endl
                << std::endl
                << " DavidsonSolver: convergence achieved..." << std::endl;
        retval = 0;
    }
    return retval;
}

// template <class OrbitalsType, class MatrixType>
// void DavidsonSolver<OrbitalsType,MatrixType>::swapColumnsVect(
//    MatrixType& evect,
//    const MatrixType& hb2N,
//    const std::vector<DISTMATDTYPE>& eval,
//    MatrixType& work2N)
//{
//    const int two_numst = (int)eval.size();
//    const int numst     = two_numst / 2;
//
//    Control& ct = *(Control::instance());
//
//    // count number of values smaller than tol
//    int n            = 0;
//    const double tol = 1.e-8;
//    while (eval[n] < tol)
//    {
//        n++;
//    }
// #ifdef USE_MPI
//    MPI_Bcast(&n, 1, MPI_INT, 0, comm_);
// #endif
//    if (n <= numst) return;
//
//    // build H matrix in basis of evect
//    MatrixType hrot("hrot", two_numst, two_numst);
//    hrot = hb2N;
//    rotateSym(hrot, evect, work2N);
//
//    // compute eigenvalues of rotated H
//    std::vector<DISTMATDTYPE> hdval(two_numst);
//    hrot.getDiagonalValues(&hdval[0]);
//
//    while (n > numst)
//    {
//        int index           = n - 1;
//        DISTMATDTYPE minval = hdval[index];
//
//        // find smallest diagonal values for hrot
//        for (int i = 0; i < n - 1; i++)
//        {
//            if (hdval[i] < minval)
//            {
//                minval = hdval[i];
//                index  = i;
//            }
//        }
// #ifdef USE_MPI
//        MPI_Bcast(&index, 1, MPI_INT, 0, comm_);
// #endif
//        if (index != n - 1)
//        {
//            if (onpe0 && ct.verbose > 2)
//                os_ << "Swap vectors " << index << " and " << n - 1
//                    << " with eigenvalues (for H) " << hdval[index] << " and "
//                    << hdval[n - 1] << std::endl;
//            evect.swapColumns(index, n - 1);
//            const DISTMATDTYPE tmp = hdval[index];
//            hdval[index]           = hdval[n - 1];
//            hdval[n - 1]           = tmp;
//        }
//
//        n--;
//    }
//}

template <class OrbitalsType, class MatrixType>
double DavidsonSolver<OrbitalsType, MatrixType>::evaluateDerivative(
    MatrixType& dm2Ninit, MatrixType& delta_dm, const double ts0)
{
    work2N_->symm('l', 'l', 1., proj_mat2N_->getMatHB(), delta_dm, 0.);

    // factor 0.5 to get value in Hartree (HB is in Rydberg)
    double de0      = 0.5 * work2N_->trace();
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (mmpi.nspin() > 1)
    {
        double tmp = 0.;
        mmpi.allreduceSpin(&de0, &tmp, 1, MPI_SUM);
        de0 = tmp;
    }

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
        os_ << "Derivative of U in 0 = " << de0 << std::endl;

    //
    // evaluate numerical derivative of entropy in beta=0
    //
    // if( onpe0 ) os_<<"evaluate numerical derivative of entropy in
    // beta=0"<<endl;
    const double dbeta = 0.0001;
    *work2N_           = dm2Ninit;
    work2N_->axpy(dbeta, delta_dm);

    proj_mat2N_->setDM(*work2N_);
    proj_mat2N_->computeOccupationsFromDM();

    const double tsd0e = evalEntropy(proj_mat2N_.get(), false, os_);
    const double dts0e = (tsd0e - ts0) / dbeta;
    if (onpe0 && ct.verbose > 2)
    {
        os_ << "Numerical derivative of -TS in 0 = " << std::scientific
            << -dts0e << std::endl;
        os_ << "interval for finite differences = " << dbeta << std::endl;
    }
    de0 -= dts0e;

    return de0;
}

template <class OrbitalsType, class MatrixType>
void DavidsonSolver<OrbitalsType, MatrixType>::buildTarget2N_MVP(
    MatrixType& h11, MatrixType& h12, MatrixType& h21, MatrixType& h22,
    MatrixType& s11, MatrixType& s22, MatrixType& target)
{
    target_tm_.start();

    Control& ct = *(Control::instance());

    proj_mat2N_->assignBlocksH(h11, h12, h21, h22);

    // if( onpe0 )os_<<"Build S2N..."<<endl;
    MatrixType s2N("s2N", 2 * numst_, 2 * numst_);
    s2N.assign(s11, 0, 0);
    // s2N.assign(s12,0,numst_);
    // s2N.assign(s21,numst_,0);
    s2N.assign(s22, numst_, numst_);

    int orbitals_index = 0;
    proj_mat2N_->setGramMatrix(s2N, orbitals_index);

    //
    // compute target density matrix, with occupations>0 in numst only
    //
    // if( onpe0 )os_<<"Compute X2N..."<<endl;
    proj_mat2N_->setHB2H();

    proj_mat2N_->updateDM();

    target = proj_mat2N_->dm();

    if (ct.verbose > 2)
    {
        double pnel = proj_mat2N_->getNel();
        if (onpe0) os_ << "Number of electrons = " << pnel << std::endl;
    }

    target_tm_.stop();
}

// template <class OrbitalsType, class MatrixType>
// void DavidsonSolver<OrbitalsType,MatrixType>::buildTarget2N_new(
//    MatrixType& h11,
//    MatrixType& h12,
//    MatrixType& h21,
//    MatrixType& h22,
//    MatrixType& s11,
//    MatrixType& s22,
//    const std::vector<DISTMATDTYPE>& occ,
//    const std::vector<DISTMATDTYPE>& auxenergies, const double kbT,
//    const double eta, MatrixType& target)
//{
//    proj_mat2N_->assignBlocksH(h11, h12, h21, h22);
//
//    // if( onpe0 )os_<<"Build S2N..."<<endl;
//    MatrixType s2N("s2N", 2 * numst_, 2 * numst_);
//    s2N.assign(s11, 0, 0);
//    // s2N.assign(s12,0,numst_);
//    // s2N.assign(s21,numst_,0);
//    s2N.assign(s22, numst_, numst_);
//
//    int orbitals_index = 0;
//    proj_mat2N_->setGramMatrix(s2N, orbitals_index);
//
//    //
//    // compute target density matrix, with occupations>0 in numst only
//    //
//    // if( onpe0 )os_<<"Compute X2N..."<<endl;
//    std::vector<DISTMATDTYPE> eigenval(2 * numst_);
//    proj_mat2N_->setHB2H();
//    proj_mat2N_->solveGenEigenProblem(*work2N_, eigenval);
//    // if( onpe0 )
//    // for(int i=0;i<2*numst;i++)
//    //    os_<<"eigenvalue of HB"<<i<<"="<<eigenval[i]<<endl;
//
//    std::vector<DISTMATDTYPE> dedf(numst_, 0.);
//    const double tol = 1.e-12;
//    double mu        = 0.;
//    double mup       = 0.;
//    for (int i = 0; i < numst_; i++)
//    {
//        os_ << std::scientific;
//        os_ << std::setprecision(10);
//        // if( onpe0 )os_<<"auxenergies[i]="<<auxenergies[i]<<endl;
//        // if( onpe0 )os_<<"occ[i]="<<occ[i]<<endl;
//        // if( onpe0 )os_<<"log(1.-occ[i])="<<log(1.-occ[i])<<endl;
//        // if( onpe0 )os_<<"log(occ[i])="<<log(occ[i])<<endl;
//        if (occ[i] > tol && occ[i] < (1. - tol))
//            dedf[i]
//                = -2. * (1. - occ[i]) * occ[i]
//                  * (auxenergies[i] - kbT * (log(1. - occ[i]) - log(occ[i])))
//                  / kbT;
//        mu += dedf[i];
//        mup += -2. * (1. - occ[i]) * occ[i] / kbT;
//    }
//    if (fabs(mup) > 1.e-6) mu /= mup;
//    if (onpe0)
//    {
//        os_ << "mu=" << mu << std::endl;
//        os_ << "mup=" << mup << std::endl;
//    }
//    for (int i = 0; i < numst_; i++)
//        dedf[i] += 2. * mu * (1. - occ[i]) * occ[i] / kbT;
//
//    std::vector<DISTMATDTYPE> occ2N(2 * numst_, 0.);
//
//    if (onpe0)
//    {
//        os_ << "eta=" << eta << std::endl;
//    }
//    double sum = 0.;
//    for (int i = 0; i < numst_; i++)
//    {
//        occ2N[numst_ - i - 1] = occ[i] + eta * dedf[i];
//        sum += dedf[i];
//        if (onpe0)
//            os_ << "occ2N[numst-i-1]=" << occ2N[numst_ - i - 1] << std::endl;
//    }
//    if (onpe0) os_ << "sum=" << sum << std::endl;
//
//    // proj_mat2N_->updateAuxilliaryEnergies();
//    //
//    proj_mat2N_->computeChemicalPotentialAndOccupations(ct.occ_width,ct.getNel(),numst_);
//    proj_mat2N_->buildDM(*work2N_, occ2N, orbitals_index);
//
//    Control& ct = *(Control::instance());
//    if (ct.verbose > 2)
//    {
//        double pnel = proj_mat2N_->getNel();
//        if (onpe0) os_ << "Number of electrons = " << pnel << std::endl;
//    }
//
//    target = proj_mat2N_->dm();
//}

// update density matrix in 2N x 2N space
template <class OrbitalsType, class MatrixType>
int DavidsonSolver<OrbitalsType, MatrixType>::solve(
    OrbitalsType& orbitals, OrbitalsType& work_orbitals)
{
    assert(numst_ == static_cast<int>(orbitals.numst()));

    solve_tm_.start();

    int retval = 1; // 0 -> converged

    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (onpe0 && ct.verbose > 1)
    {
        os_ << "---------------------------------------------------------------"
               "-"
            << std::endl;
        os_ << "Update DM and wave functions using Davidson Solver..."
            << std::endl;
        os_ << "---------------------------------------------------------------"
               "-"
            << std::endl;
    }
    // step for numerical derivative

    Potentials& pot = hamiltonian_->potential();

    // double nel=orbitals.projMatrices()->getNel();
    // if( onpe0 )
    //    os_<<"NEW algorithm: Number of electrons at start = "<<nel<<endl;
    de_old_ = 100000.;
    de_     = 100000.;

    KBPsiMatrixSparse kbpsi_1(nullptr);
    kbpsi_1.setup(ions_);
    KBPsiMatrixSparse kbpsi_2(nullptr);
    kbpsi_2.setup(ions_);

    MatrixType s11("s11", numst_, numst_);
    MatrixType s22("s22", numst_, numst_);
    s11.identity();
    s22.identity();

    // orthonormalize initial wave functions to avoid drifting away from
    // orthonormalization over outer iterations (MD or geometry optimization)
    // jlf, 01/02/2021
    if (mmpi.PE0() && ct.verbose > 1)
        os_ << "Orthonormalize wavefunctions at start of DavidsonSolver"
            << std::endl;
    orbitals.orthonormalizeLoewdin(false, nullptr, false);

    for (int outer_it = 0; outer_it <= ct.max_electronic_steps; outer_it++)
    {
        // turn on PB solver if necessary
        diel_control_.activate(outer_it, de_);

        ProjectedMatrices<MatrixType>* proj_matN
            = dynamic_cast<ProjectedMatrices<MatrixType>*>(
                orbitals.getProjMatrices());
        if (mmpi.PE0() && ct.verbose > 0)
        {
            os_ << "###########################" << std::endl;
            os_ << "DavidsonSolver -> Iteration " << outer_it << std::endl;
            os_ << "###########################" << std::endl;
        }
        OrbitalsType hphi("Davidson_hphi", orbitals);
        MatrixType dm2Ninit("dm2N", 2 * numst_, 2 * numst_);
        std::vector<DISTMATDTYPE> eval(2 * numst_);
        MatrixType evect("EigVect", 2 * numst_, 2 * numst_);

        MatrixType dm11("dm11", numst_, numst_);
        MatrixType dm12("dm12", numst_, numst_);
        MatrixType dm21("dm21", numst_, numst_);
        MatrixType dm22("dm22", numst_, numst_);

        // save computed vh for a fair energy "comparison" with vh computed
        // in close neigborhood
        const pb::GridFunc<POTDTYPE> vh_init(electrostat_->getVh());

        orbitals.setDataWithGhosts();

        // Update density
        rho_->update(orbitals);

        MatrixType h11("h11", numst_, numst_);
        MatrixType h22("h22", numst_, numst_);
        MatrixType h12("h12", numst_, numst_);
        MatrixType h21("h21", numst_, numst_);

        MatrixType h11nl("h11nl", numst_, numst_);
        MatrixType h22nl("h22nl", numst_, numst_);
        MatrixType h12nl("h12nl", numst_, numst_);

        kbpsi_1.computeAll(ions_, orbitals);

        kbpsi_1.computeHvnlMatrix(&kbpsi_1, ions_, h11nl);

        for (int inner_it = 0; inner_it < ct.dm_inner_steps; inner_it++)
        {
            if (mmpi.PE0() && ct.verbose > 1)
            {
                os_ << "---------------------------" << std::endl;
                os_ << "Inner iteration " << inner_it << std::endl;
                os_ << "---------------------------" << std::endl;
            }

            // Update potential
            mgmol_strategy_->update_pot(vh_init, ions_);

            energy_->saveVofRho();

            ProjectedMatrices<MatrixType>* current_proj_mat
                = (inner_it == 0) ? proj_matN : proj_mat2N_.get();
            if (ct.verbose > 2) current_proj_mat->printOccupations(os_);

            double ts0       = 0.;
            double e0        = 0.;
            const int printE = (ct.verbose > 1 || outer_it % 10 == 0) ? 1 : 0;
            if (inner_it == 0)
            {
                // orbitals are new, so a few things need to recomputed
                ProjectedMatrices<MatrixType>* projmatrices
                    = dynamic_cast<ProjectedMatrices<MatrixType>*>(
                        orbitals.getProjMatrices());
                assert(projmatrices != nullptr);

                // get H*phi stored in hphi
                // h11 computed at the same time
                mgmol_strategy_->computePrecondResidual(orbitals, hphi,
                    work_orbitals, ions_, &kbpsi_1, false, false);

                projmatrices->setHB2H();

                h11 = projmatrices->getH();

                // if( onpe0 )os_<<"H11..."<<endl;
                // h11.print(os_,0,0,5,5);
                // if( onpe0 )os_<<"Matrices..."<<endl;
                // projmatrices->printMatrices(os_);

                ts0 = evalEntropy(projmatrices, (ct.verbose > 1), os_);
                e0  = energy_->evaluateTotal(
                    ts0, projmatrices, ions_, orbitals, printE, os_);

                retval = checkConvergence(e0, outer_it, ct.conv_tol);
                if (retval == 0 || (outer_it == ct.max_electronic_steps))
                {
                    break;
                }

                orbitals.projectOut(work_orbitals);
                // normalize first as these vectors can become quite small...
                work_orbitals.normalize();
                work_orbitals.orthonormalizeLoewdin(false, nullptr, false);
                work_orbitals.setIterativeIndex(
                    orbitals.getIterativeIndex() + 1000);

                work_orbitals.setDataWithGhosts();
                kbpsi_2.computeAll(ions_, work_orbitals);

                kbpsi_2.computeHvnlMatrix(&kbpsi_2, ions_, h22nl);
                kbpsi_1.computeHvnlMatrix(&kbpsi_2, ions_, h12nl);

                h12 = h12nl;
                h22 = h22nl;
            }
            else
            {
                hamiltonian_->applyDeltaPot(orbitals, hphi);
                orbitals.addDotWithNcol2Matrix(hphi, h11);
            }

            if (inner_it == 0)
            {
                // compute H*P and store in hphi
                hamiltonian_->applyLocal(numst_, work_orbitals, hphi);
            }
            else
            {
                hamiltonian_->applyDeltaPot(work_orbitals, hphi);
            }

            // update h22, h12 and h21
            orbitals.addDotWithNcol2Matrix(hphi, h12);

            work_orbitals.addDotWithNcol2Matrix(hphi, h22);

            h21.transpose(1., h12, 0.);

            if (inner_it == 0)
            {
                dm2Ninit.assign(current_proj_mat->dm(), 0, 0);
            }
            else
            {
                dm2Ninit = proj_mat2N_->dm();

                proj_mat2N_->assignBlocksH(h11, h12, h21, h22);
                proj_mat2N_->setHB2H();

                ts0 = evalEntropy(proj_mat2N_.get(), (ct.verbose > 1), os_);
                e0  = energy_->evaluateTotal(
                    ts0, proj_mat2N_.get(), ions_, orbitals, printE, os_);
            }

            // 2N x 2N target...
            proj_mat2N_->setHiterativeIndex(
                orbitals.getIterativeIndex(), pot.getIterativeIndex());

            MatrixType target("target", 2 * numst_, 2 * numst_);

            MatrixType delta_dm("delta_dm", 2 * numst_, 2 * numst_);

            buildTarget2N_MVP(h11, h12, h21, h22, s11, s22, target);
            if (mmpi.instancePE0() && ct.verbose > 1 && (outer_it % 10 == 0))
            {
                proj_mat2N_->printEigenvalues(os_);
                proj_mat2N_->printOccupations(os_);
            }
            delta_dm = target;
            delta_dm -= dm2Ninit;
            double beta = 0.;
            if (mixing_ > 0.)
            {
                beta = mixing_;
                if (mmpi.PE0() && ct.verbose > 1)
                    os_ << "Davidson with beta = " << beta << std::endl;
            }
            else
            {
                double de0 = evaluateDerivative(dm2Ninit, delta_dm, ts0);

                //
                // evaluate free energy at beta=1
                //
                if (mmpi.PE0() && ct.verbose > 2)
                    os_ << "Target energy..." << std::endl;
                proj_mat2N_->setDM(target);
                proj_mat2N_->computeOccupationsFromDM();
                double nel = proj_mat2N_->getNel();
                if (mmpi.PE0() && ct.verbose > 2)
                    os_ << "Number of electrons at beta=1 : " << nel
                        << std::endl;

                dm11.getsub(target, numst_, numst_, 0, 0);
                dm12.getsub(target, numst_, numst_, 0, numst_);
                dm21.getsub(target, numst_, numst_, numst_, 0);
                dm22.getsub(target, numst_, numst_, numst_, numst_);

                rho_->computeRho(
                    orbitals, work_orbitals, dm11, dm12, dm21, dm22);

                mgmol_strategy_->update_pot(vh_init, ions_);

                energy_->saveVofRho();

                // update h11, h22, h12, and h21
                hamiltonian_->applyDeltaPot(orbitals, hphi);
                orbitals.addDotWithNcol2Matrix(hphi, h11);

                hamiltonian_->applyDeltaPot(work_orbitals, hphi);
                work_orbitals.addDotWithNcol2Matrix(hphi, h22);
                orbitals.addDotWithNcol2Matrix(hphi, h12);

                h21.transpose(1., h12, 0.);

                // assemble 2N x 2N Hamiltonian
                proj_mat2N_->assignBlocksH(h11, h12, h21, h22);
                proj_mat2N_->setHB2H();

                const double ts1
                    = evalEntropy(proj_mat2N_.get(), (ct.verbose > 2), os_);
                const double e1 = energy_->evaluateTotal(ts1, proj_mat2N_.get(),
                    ions_, orbitals, ct.verbose - 1, os_);

                // line minimization
                beta = minQuadPolynomial(e0, e1, de0, (ct.verbose > 2), os_);

                if (mmpi.PE0() && ct.verbose > 0)
                {
                    os_ << std::setprecision(12);
                    os_ << "ts1=" << ts1 << std::endl;
                    os_ << std::fixed << "Inner iteration " << inner_it
                        << ", E0=" << e0 << ", E1=" << e1;
                    os_ << std::scientific << ", E0'=" << de0
                        << " -> beta=" << beta;
                    // os_<<scientific<<", E1'="<<de1;
                    os_ << std::endl;
                }
            }

            // update DM
            *work2N_ = dm2Ninit;
            work2N_->axpy(beta, delta_dm);
            proj_mat2N_->setDM(*work2N_);

            if (inner_it < ct.dm_inner_steps - 1)
            {
                dm11.getsub(*work2N_, numst_, numst_, 0, 0);
                dm12.getsub(*work2N_, numst_, numst_, 0, numst_);
                dm21.getsub(*work2N_, numst_, numst_, numst_, 0);
                dm22.getsub(*work2N_, numst_, numst_, numst_, numst_);

                if (ct.verbose > 2)
                {
                    double pnel = proj_mat2N_->getNel();
                    if (mmpi.PE0())
                        os_ << "Number of electrons for interpolated DM = "
                            << pnel << std::endl;
                }

                rho_->computeRho(
                    orbitals, work_orbitals, dm11, dm12, dm21, dm22);
            }

        } // inner iterations

        // update orbitals
        proj_mat2N_->diagonalizeDM(eval, evect);
#if 1
        if (mmpi.PE0() && ct.verbose > 2)
        {
            os_ << "Eigenvalues of Interpolated DM: " << std::endl;
            os_ << std::fixed << std::setprecision(4);
            for (unsigned int i = 0; i < eval.size(); i++)
            {
                os_ << eval[i];
                if (i % 10 == 9)
                    os_ << std::endl;
                else
                    os_ << "  ";
            }
            os_ << std::endl;
            double tot = 0.;
            for (int i = 0; i < numst_; i++)
            {
                tot += eval[numst_ + i];
            }
            os_ << "Total occupations for top half states="
                << std::setprecision(15) << tot << std::endl;
        }
        if (mmpi.PE0() && ct.verbose > 1)
        {
            os_ << std::setprecision(15)
                << "Last level occupancy = " << eval[numst_] << std::endl;
        }
#endif

        // to differentiate very low occupation vectors, extract those with
        // lowest energy
        //        MatrixType  z2N("z2N",2*numst_,
        //        2*numst_); swapColumnsVect(evect, proj_mat2N_->getMatHB(),
        //        eval, z2N );

        dm12.getsub(evect, numst_, numst_, 0, numst_);
        dm22.getsub(evect, numst_, numst_, numst_, numst_);

        // replace orbitals with eigenvectors corresponding to largest
        // eigenvalues of DM
        orbitals.multiply_by_matrix(dm12);
        work_orbitals.multiply_by_matrix(dm22);
        orbitals.axpy((ORBDTYPE)1., work_orbitals);
        orbitals.incrementIterativeIndex();
        orbitals.incrementIterativeIndex();
        work_orbitals.incrementIterativeIndex(2);

        std::vector<double> new_occ(numst_);
        double tocc              = 0.;
        const double spin_factor = mmpi.nspin() > 1 ? 1. : 0.5;
        for (int i = 0; i < numst_; i++)
        {
            new_occ[i] = spin_factor * eval[numst_ + i];
            tocc += new_occ[i];
        }
        double kbT = ct.occ_width;
        if (mmpi.PE0() && ct.verbose > 2)
            os_ << "Total occupations/spin for kbT = " << kbT << ": "
                << std::setprecision(15) << tocc << std::endl;

#if 0 // use lower T?
        const double tol_sum_occ=1.e-2;
        if( (0.5*ct.getNel()-tocc)>tol_sum_occ )
        {
            ProjectedMatrices* projmatrices=dynamic_cast<ProjectedMatrices*>( orbitals.getProjMatrices() );
            assert( projmatrices!=0 );
            projmatrices->setOccupations(new_occ);
            projmatrices->setAuxilliaryEnergiesFromOccupations();
            while( (0.5*ct.getNel()-tocc)>tol_sum_occ )
            {
                kbT*=0.5;
                projmatrices->computeChemicalPotentialAndOccupations(kbT,ct.getNel(),numst_);
                projmatrices->getOccupations(new_occ);
            
                tocc=0.;
                for(int j=0;j<numst_;j++)
                {
                    tocc+=new_occ[j];
                }
                if( onpe0 && ct.verbose>2 )
                    os_<<"Total occupations for kbT = "<<kbT<<": "<<setprecision(15)<<2.*tocc<<endl;
            }
        }
        
        //if( onpe0 && (it%10)==0 && ct.verbose>2 )
        if( onpe0 )
        {
            os_<<"New occupations: "<<endl;
            os_<<fixed<<setprecision(4);
            for(int i=0;i<numst_;i++)
            {
                os_<<new_occ[i];
                if( i%10==9 )os_<<endl;
                else         os_<<"  ";
            }
            os_<<endl;
        }
#endif
#if 1
        // add occupation to lowest states until correct number of e- is reached
        if (mmpi.PE0() && ct.verbose > 2)
            os_ << "Total occupations/spin before correction = "
                << std::setprecision(15) << tocc << std::endl;
        int j                   = numst_ - 1;
        const double target_nel = ct.getNelSpin();
        while ((target_nel - tocc) > 1.e-8 && j >= 0)
        {
            const double delta = target_nel - tocc;
            tocc -= new_occ[j];
            new_occ[j] = std::min(1., new_occ[j] + delta);
            tocc += new_occ[j];
            j--;
            // if( onpe0 && ct.verbose>2 )
            //    os_<<"Total occupations = "<<setprecision(15)<<2.*tocc<<endl;
        }
        if (mmpi.PE0() && ct.verbose > 2)
            os_ << "Total occupations/spin = " << std::setprecision(15) << tocc
                << std::endl;

        if (mmpi.PE0() && ((outer_it % 10) == 0 || (retval == 0))
            && ct.verbose > 2)
        {
            os_ << "Final Occupations: " << std::endl;
            os_ << std::fixed << std::setprecision(4);
            for (int i = 0; i < numst_; i++)
            {
                os_ << new_occ[i];
                if (i % 10 == 9)
                    os_ << std::endl;
                else
                    os_ << "  ";
            }
            os_ << std::endl;
        }
#endif

        // build diagonal DM (orbitals are now eigenvectors of DM)
        ProjectedMatrices<MatrixType>* pmat
            = dynamic_cast<ProjectedMatrices<MatrixType>*>(
                orbitals.getProjMatrices());
        assert(pmat);
        pmat->buildDM(new_occ);

        if (retval == 0) break;

    } // loop over wave functions updates

    // Generate new density
    rho_->update(orbitals);

    if (mmpi.PE0() && ct.verbose > 0)
    {
        ProjectedMatrices<MatrixType>* pmat
            = dynamic_cast<ProjectedMatrices<MatrixType>*>(
                orbitals.getProjMatrices());
        assert(pmat);

        pmat->printOccupations(os_);
        proj_mat2N_->printEigenvalues(os_);
    }

    if (mmpi.PE0() && ct.verbose > 1)
    {
        os_ << "------------------------------------------------------------"
            << std::endl;
        os_ << "End Davidson Solver..." << std::endl;
        os_ << "------------------------------------------------------------"
            << std::endl;
    }
    solve_tm_.stop();

    return retval;
}

template <class OrbitalsType, class MatrixType>
void DavidsonSolver<OrbitalsType, MatrixType>::printTimers(std::ostream& os)
{
    solve_tm_.print(os);
    target_tm_.print(os);
}

#ifdef MGMOL_USE_SCALAPACK
template class DavidsonSolver<ExtendedGridOrbitals<ORBDTYPE>,
    dist_matrix::DistMatrix<DISTMATDTYPE>>;
#endif
template class DavidsonSolver<ExtendedGridOrbitals<ORBDTYPE>, ReplicatedMatrix>;
