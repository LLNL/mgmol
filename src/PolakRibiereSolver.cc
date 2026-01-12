// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "PolakRibiereSolver.h"

#include "Control.h"
#include "DMStrategy.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "Potentials.h"
#include "ProjectedMatricesInterface.h"
#include "Rho.h"

template <class OrbitalsType>
Timer PolakRibiereSolver<OrbitalsType>::solve_tm_("solve");
template <class OrbitalsType>
int PolakRibiereSolver<OrbitalsType>::it_scf_ = 0;

template <class OrbitalsType>
PolakRibiereSolver<OrbitalsType>::PolakRibiereSolver(
    Hamiltonian<OrbitalsType>* hamiltonian,
    ProjectedMatricesInterface* proj_matrices, Energy<OrbitalsType>* energy,
    Electrostatic* electrostat, MGmol<OrbitalsType>* mgmol_strategy, Ions& ions,
    Rho<OrbitalsType>* rho, DMStrategy<OrbitalsType>* dm_strategy,
    std::ostream& os)
    : mgmol_strategy_(mgmol_strategy),
      hamiltonian_(hamiltonian),
      proj_matrices_(proj_matrices),
      energy_(energy),
      electrostat_(electrostat),
      ions_(ions),
      rho_(rho),
      dm_strategy_(dm_strategy),
      os_(os)
{
    Control& ct(*(Control::instance()));

    alpha_ = ct.precond_factor;
    r_k_   = nullptr;
    r_km1_ = nullptr;

    p_k_ = nullptr;

    z_k_   = nullptr;
    z_km1_ = nullptr;

    sigma_a_ = 1.e-4;
    sigma_b_ = 1.e-1;

    with_preconditioner_ = (ct.withPreconditioner());
}

template <class OrbitalsType>
PolakRibiereSolver<OrbitalsType>::~PolakRibiereSolver()
{
}

template <class OrbitalsType>
void PolakRibiereSolver<OrbitalsType>::printEnergy(const short step) const
{
    if (onpe0)
    {
        os_ << std::setprecision(12) << std::fixed << "%%    " << step
            << " SC ENERGY = " << eks_history_[0];
        if (step > 0)
        {
            if (testUpdatePot())
            {
                os_ << std::setprecision(2) << std::scientific
                    << ", delta Eks = " << de_ << std::endl;
            }
            else
            {
                os_ << std::setprecision(2) << std::scientific
                    << ", delta Eig. sum = " << deig_ << std::endl;
                os_ << std::setprecision(12) << std::fixed
                    << "Sum eigenvalues = " << sum_eig_[0] << std::endl;
            }
        }
        else
        {
            os_ << std::endl;
        }
    }
}

template <class OrbitalsType>
bool PolakRibiereSolver<OrbitalsType>::checkPrintResidual(
    const short step) const
{
    Control& ct(*(Control::instance()));
    return (ct.iprint_residual > 0) ? !(step % ct.iprint_residual) : false;
}

template <class OrbitalsType>
void PolakRibiereSolver<OrbitalsType>::dielON()
{
    Potentials& pot(hamiltonian_->potential());
    bool isON = pot.diel();
    if (!isON) return; // continuum solvent is OFF

    static bool pbset = false;
    if (pbset) return; // continuum solvent already set

    Control& ct(*(Control::instance()));
    const int diel_delay = (ct.restart_info < 3) ? 10 : 1;
    if (it_scf_ < diel_delay)
    {
        isON = false;
    }

    // turn ON continuum solvent
    if (isON && (ct.restart_info >= 3 || deig2_ < 2.e-2 * ct.numst))
    {
        electrostat_->setupPB(ct.e0_, ct.rho0_, ct.drho0_, pot);
        electrostat_->setup(ct.vh_its);
        pbset = true;
    }

    if (pot.diel() && !pbset)
        if (onpe0 && ct.verbose > 1)
            os_ << " Solvation turned off for this step" << std::endl;
}

template <class OrbitalsType>
bool PolakRibiereSolver<OrbitalsType>::testUpdatePot() const
{
    Control& ct(*(Control::instance()));
    return (it_scf_ > ct.max_changes_pot);
}

template <class OrbitalsType>
bool PolakRibiereSolver<OrbitalsType>::checkConvPot() const
{
    Control& ct(*(Control::instance()));
    Potentials& pot(hamiltonian_->potential());
    if ((fabs(pot.scf_dvrho() / ions_.getNumIons()) < ct.conv_tol)
        && (it_scf_ > (2 + ct.max_changes_pot)))
        return true;

    return false;
}

template <class OrbitalsType>
bool PolakRibiereSolver<OrbitalsType>::checkWolfeConditions(
    const double trial_step_energy, const double alpha_k) const
{
    assert(sigma_a_ > 0.);
    assert(sigma_b_ > 0.);

    // const double dk  =-1.*r_k_->dotProduct(*p_k_);
    const double dkm1 = -1.
                        * r_km1_->dotProduct(*p_k_,
                            2); // inverse(S) already included in r_km1_

    const bool wolfe0
        = (trial_step_energy <= eks_history_[0] + sigma_a_ * alpha_k * dkm1);
    // const bool wolfe1 = ( dk >= sigma_b_*dkm1 );
    const bool wolfe1 = true;

    os_ << std::setprecision(10);
    if (!wolfe0 && onpe0)
    {
        os_ << "wolfe0 is not satisfied:" << std::endl;
        os_ << "trial_step_energy                     =" << trial_step_energy
            << std::endl;
        os_ << "eks_history_[0]=" << eks_history_[0] << std::endl;
        os_ << "eks_history_[0]+sigma_a_*alpha_k*dkm1="
            << eks_history_[0] + sigma_a_ * alpha_k * dkm1 << std::endl;
    }
    // if( !wolfe1 && onpe0 )os_<<"wolfe1 is not satisfied: dk="<<dk
    //                         <<", dkm1="<<dkm1<<endl;

    return wolfe0 && wolfe1;
}

// returns:
// 0 if converged,
// 1 if not converged yet,
// -1 if not reaching minimum convergence
// -2 if failing to converge
template <class OrbitalsType>
int PolakRibiereSolver<OrbitalsType>::checkConvergenceEnergy(
    const short step, const short max_steps)
{
    Control& ct(*(Control::instance()));

    double dtol2 = 1000.;
    if (step > 0)
    {
        double deig_old = deig_;
        deig_           = std::abs(sum_eig_[1] - sum_eig_[0]);
        deig2_          = std::max(deig_, deig_old);

        // energy variation during last iteration
        double de_old = de_;
        de_           = std::abs(eks_history_[0] - eks_history_[1]);
        de2_          = std::max(de_, de_old);

        if (testUpdatePot())
            dtol2 = de2_;
        else
            dtol2 = deig2_;
    }

    // test if this iteration converged
    if (step > 1 && dtol2 < ct.conv_tol) return 0;

    if (step == max_steps && dtol2 > ct.conv_tol_stop) return -1;

    // Test very bad behavior!!
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    Potentials& pot(hamiltonian_->potential());
    if (pot.scf_dv() > 1. && step > 50 && step > ct.max_changes_pot)
    {
        if (onpe0)
        {
            os_ << std::endl
                << "MGmol ERROR, electronic iterations are not converging at "
                   "all: DV ="
                << pot.scf_dv() << std::endl
                << std::endl;
            os_ << std::flush;
        }
        mmpi.barrier();
        return -2;
    }

    return 1;
}

template <class OrbitalsType>
double PolakRibiereSolver<OrbitalsType>::evaluateEnergy(
    const OrbitalsType& orbitals, const bool print_flag)
{
    // Get the new total energy
    const double ts     = 0.5 * proj_matrices_->computeEntropy(); // in [Ha]
    const double energy = energy_->evaluateTotal(
        ts, proj_matrices_, ions_, orbitals, print_flag, os_);

    // Control& ct(*(Control::instance()));
    // if( ct.verbose>2 && onpe0 )os_<<"energy="<<energy<<endl;

    sum_eig_[1] = sum_eig_[0];
    sum_eig_[0] = 2. * proj_matrices_->getEigSum(); // 2.*sum in [Ry]

    return energy;
}

// Polak-Ribiere
template <class OrbitalsType>
double PolakRibiereSolver<OrbitalsType>::computeBeta(
    OrbitalsType& work_orbitals) const
{
    work_orbitals.assign(*r_k_);
    work_orbitals.axpy((ORBDTYPE)(-1.), *r_km1_);

    double beta = z_k_->dotProduct(work_orbitals, 2);

    beta /= z_km1_->dotProduct(*r_km1_, 2);

    return beta;
}

template <class OrbitalsType>
int PolakRibiereSolver<OrbitalsType>::solve(OrbitalsType& orbitals,
    OrbitalsType& work_orbitals, Ions& ions, const short max_steps,
    const short iprint, double& last_eks)
{
    Control& ct(*(Control::instance()));
    assert(ct.getOrthoType() != OrthoType::Eigenfunctions);

    solve_tm_.start();

    eks_history_[0] = 1.e9;
    eks_history_[1] = 1.e9;

    int retval = 1; // 0 -> converged, -1 -> problem, -2 -> ( de>conv_tol_stop )

    deig_  = 1.e9;
    deig2_ = 1.e9;
    de_    = 1.e9;
    de2_   = 1.e9;

#ifdef HAVE_ARPACK
    if (ct.precond_factor_computed)
    {
        const double ela = getLAeigen(0., 500, ions);
        if (fabs(ela) > 1.e-16)
            ct.precond_factor = 1. / ela;
        else
            return -1;
    }
#else
    if (ct.precond_factor_computed)
    {
        os_ << "Needs ARPACK to compute Preconditioner factor" << std::endl;
        return -1;
    }
#endif
    if (onpe0)
    {
        os_ << "### PolakRibiereSolver Polak-Ribiere###" << std::endl;
        if (ct.verbose > 1 && ct.precond_factor_computed)
            os_ << "Preconditioning factor: " << ct.precond_factor << std::endl;
    }

    r_k_   = new OrbitalsType("PR_rk", orbitals, false);
    r_km1_ = new OrbitalsType("PR_rkm1", orbitals, false);

    if (with_preconditioner_)
    {
        z_k_   = new OrbitalsType("PR_zk", orbitals, false);
        z_km1_ = new OrbitalsType("PR_zkm1", orbitals, false);
    }
    p_k_ = new OrbitalsType("PR_pk", orbitals, false);

    int dm_success       = 0;
    bool wolfe           = false;
    bool skip_wolfe_cond = true;
    double alpha_k       = alpha_;

    // main electronic structure outer loop
    for (short step = 0; step <= max_steps; step++)
    {
        // if( onpe0 )os_<<"Polak-Ribiere step "<<step<<endl;

        if (dm_success == 0) // if DM determined at previous step is good
        {
            // turn on PB solver if necessary
            dielON();

            proj_matrices_->resetDotProductMatrices();

            // Generate new density
            rho_->update(orbitals);

            // Update potential
            if (testUpdatePot()) mgmol_strategy_->update_pot(ions);

            mgmol_strategy_->updateHmatrix(orbitals, ions);

            // theta = invB * Hij
            // (to be used for energy and gradient computation)
            proj_matrices_->updateThetaAndHB();

            if (step == max_steps)
            {
                solve_tm_.stop();
                break;
            }
            // Output the eigenvalues and occupations
            bool flag = false;
            if (iprint) flag = (!(step % iprint) && step < max_steps - 1);

            if (flag && !ct.short_sighted) mgmol_strategy_->printEigAndOcc();

            last_eks = evaluateEnergy(orbitals, flag);

            // bool print_res = checkPrintResidual(step);
            bool print_res = true;

            // evaluate residuals, preconditioned residuals for current orbitals
            double normRes = mgmol_strategy_->computeResidual(orbitals,
                work_orbitals, ions_, *r_k_, (print_res || ct.checkResidual()),
                ct.checkResidual());
            if (normRes < ct.conv_tol && ct.checkResidual())
            {
                solve_tm_.stop();
                return 0;
            }

            z_k_->assign(*r_k_);
            if (with_preconditioner_)
            {
                // PRECONDITIONING
                // compute the preconditioned steepest descent direction
                // -> res
                mgmol_strategy_->precond_mg(*z_k_);
            }

            // save actual residual*inverse(S) (taking into account occupations)
            // SquareLocalMatrices<MATDTYPE,MemorySpace::Host>& localX(
            // proj_matrices_->getLocalX() ); r_k_->multiplyByMatrix(localX);

            // non-energy compatible spread penalties are added after
            // preconditioning is applied
            if (ct.isSpreadFunctionalActive() && !ct.isSpreadFunctionalEnergy())
                mgmol_strategy_->addResidualSpreadPenalty(orbitals, *z_k_);

            // apply AOMM projector to have a gradient orthogonal to kernel
            // functions
            if (ct.use_kernel_functions)
            {
                mgmol_strategy_->projectOutKernel(*z_k_);

                // just reapply mask for now...
                z_k_->applyMask();
            }

            wolfe = (skip_wolfe_cond) ? true
                                      : checkWolfeConditions(last_eks, alpha_k);
        }

        // strip dm from the overlap contribution
        // dm <- Ls**T * dm * Ls
        dm_strategy_->stripDM();

        //
        // trial step to update orbitals using r_k_, z_k_
        //
        if (wolfe && dm_success == 0) // last step was accepted
        {
            // save energy recent history
            eks_history_[1] = eks_history_[0]; // save to check convergence
            eks_history_[0] = last_eks;

            // test for convergence
            retval = checkConvergenceEnergy(step, max_steps);

            // terminate if convergence problem
            if (retval < 0) return retval;

            printEnergy(step);

            // terminate early if convergence achieved
            if (retval == 0 && !ct.checkResidual())
            {
                if (onpe0)
                    os_ << std::endl
                        << std::endl
                        << " PolakRibiereSolver: convergence achieved for "
                           "delta E..."
                        << std::endl;
                solve_tm_.stop();
                break;
            }

            // update p_k_
            if (step == 0)
                p_k_->assign(*z_k_);
            else
            {
                double beta = computeBeta(work_orbitals);
                if (beta < 0.) beta = 0.;
                if (onpe0 && ct.verbose > 0)
                    os_ << " PolakRibiereSolver: beta=" << beta << std::endl;
                p_k_->scal(beta);
                p_k_->axpy((ORBDTYPE)1., *z_k_);

                if (beta > 0.1) p_k_->scal(1. / (1. + beta));
            }

            // reset alpha to original value
            alpha_k = alpha_;

            // make new trial step
            orbitals.axpy((ORBDTYPE)alpha_k, *p_k_);

            // save current "k" vectors into "km1" vectors
            if (with_preconditioner_)
                z_km1_->assign(*z_k_);
            else
                z_km1_ = r_km1_;
            r_km1_->assign(*r_k_);
        }
        else // !wolfe
        {
            alpha_k *= 0.5;
            if (onpe0)
                os_ << "PolakRibiereSolver: Wolfe conditions not satisfied, "
                       "reduce alpha to "
                    << alpha_k << "..." << std::endl;

            // half step back
            orbitals.axpy((ORBDTYPE)(-1. * alpha_k), *p_k_);

            // return DM to value before trial step
            dm_strategy_->reset();
        }

        orbitals.incrementIterativeIndex();

        if (ct.getOrthoType() != OrthoType::Orthonormal)
        {
            orbitals.normalize();
        }

        // recompute overlap
        if (ct.getOrthoType() != OrthoType::Orthonormal)
        {
            orbitals.computeGramAndInvS();
        }

#ifdef MGMOL_USE_SCALAPACK
        // rotate pairs if smallest eigenvalue of overlap matrix below threshold
        if (ct.getThresholdEigenvalueGramQuench() > 0. && wolfe)
        {
            bool wannier_pairs = mgmol_strategy_->rotateStatesPairsOverlap(
                orbitals, work_orbitals, ct.getThresholdEigenvalueGramQuench());

            if (wannier_pairs)
            {
                orbitals.computeGramAndInvS();
                // reset reference energy since rotation may be too much of a
                // perturbation
                eks_history_[0] = 1.e8;
                eks_history_[1] = 1.e8;
                if (onpe0)
                    os_ << "PolakRibiereSolver: reset reference energy"
                        << std::endl;
            }
        }
#endif

        // rebuild dm with new overlap matrix
        dm_strategy_->dressDM();

        // update ghost data now, so that it's ready to use later
        orbitals.setDataWithGhosts();
        orbitals.trade_boundaries();

        // compute B and its inverse for Mehrstellen, so that it's ready to use
        // later (depends on Phi only)
        orbitals.computeBAndInvB(*(hamiltonian_->lapOper()));

        if (ct.computeCondGramQuench())
        {
            double condS = proj_matrices_->computeCond();
            if (onpe0)
                os_ << std::setprecision(2) << std::scientific
                    << "Condition Number of S: " << condS << std::endl;
        }

        if (testUpdatePot())
        {
            // updated Hij needed to compute new DM
            if (dm_strategy_->needH())
                mgmol_strategy_->updateHmatrix(orbitals, ions);

            // compute new density matrix
            dm_success = dm_strategy_->update(orbitals);

            skip_wolfe_cond = false;
        }

        incInnerIt();

    } // end iterations

    if (r_k_ != nullptr) delete r_k_;
    if (r_km1_ != nullptr) delete r_km1_;
    if (p_k_ != nullptr) delete p_k_;
    if (z_k_ != nullptr) delete z_k_;
    if (z_km1_ != nullptr) delete z_km1_;

    if (iprint && !ct.short_sighted) mgmol_strategy_->printEigAndOcc();

#ifdef HAVE_ARPACK
    if (ct.precond_factor_computed)
    {
        const double alpha = 1. / getLAeigen(0., 500, ions);
        if (onpe0)
        {
            os_ << "Quench electrons" << endl;
            os_ << "A Posteriori Preconditioning factor: " << alpha << endl;
        }
    }
#endif
    solve_tm_.stop();

    return retval;
}

template class PolakRibiereSolver<LocGridOrbitals<ORBDTYPE>>;
template class PolakRibiereSolver<ExtendedGridOrbitals<ORBDTYPE>>;
