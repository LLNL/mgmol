// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DFTsolver.h"
#include "ABPG.h"
#include "Control.h"
#include "DMStrategy.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "Ions.h"
#include "MGmol.h"
#include "Potentials.h"
#include "ProjectedMatricesInterface.h"
#include "Rho.h"

#ifdef MGMOL_USE_SCALAPACK
#include "GrassmanCGFactory.h"
#endif

template <class OrbitalsType>
DFTsolver<OrbitalsType>::DFTsolver(Hamiltonian<OrbitalsType>* hamiltonian,
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
      os_(os),
      diel_control_(hamiltonian->potential(), electrostat, os, onpe0)
{
    Control& ct(*(Control::instance()));

    switch (ct.OuterSolver())
    {
        case OuterSolverType::ABPG:
        {
            orbitals_stepper_
                = new ABPG<OrbitalsType>(hamiltonian_, mgmol_strategy, os_);
            break;
        }

#ifdef MGMOL_USE_SCALAPACK
        case OuterSolverType::NLCG:
        {
            orbitals_stepper_ = GrassmanCGFactory<OrbitalsType>::create(
                hamiltonian_, proj_matrices_, mgmol_strategy, ions, os_,
                ct.short_sighted);

            break;
        }
#endif

        default:
            std::cerr << "DFTsolver: Undefined iterative electronic structure "
                         "solver!!!"
                      << std::endl;
    }

    accelerate_ = false;
    // if(ct.restart_info>2 && ct.atoms_dyn == 2 && ct.wf_dyn==1
    // )accelerate_=true;
}

template <class OrbitalsType>
DFTsolver<OrbitalsType>::~DFTsolver()
{
    delete orbitals_stepper_;
}

template <class OrbitalsType>
void DFTsolver<OrbitalsType>::printEnergy(const short step) const
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
bool DFTsolver<OrbitalsType>::checkPrintResidual(const short step) const
{
    Control& ct(*(Control::instance()));
    return (ct.iprint_residual > 0) ? !(step % ct.iprint_residual) : false;
}

template <class OrbitalsType>
bool DFTsolver<OrbitalsType>::testUpdatePot() const
{
    Control& ct(*(Control::instance()));
    return (it_scf_ > ct.max_changes_pot);
}

template <class OrbitalsType>
bool DFTsolver<OrbitalsType>::checkConvPot() const
{
    Control& ct(*(Control::instance()));
    Potentials& pot(hamiltonian_->potential());
    if ((fabs(pot.scf_dvrho() / ions_.getNumIons()) < ct.conv_tol)
        && (it_scf_ > (2 + ct.max_changes_pot)))
        return true;

    return false;
}

// returns:
// 0 if converged,
// 1 if not converged yet,
// -1 if not reaching minimum convergence
// -2 if failing to converge
template <class OrbitalsType>
int DFTsolver<OrbitalsType>::checkConvergenceEnergy(
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

    // accelerate only if close enough to solution
    accelerate_ = (accelerate_ || (ct.wf_dyn && (deig2_ < 1.e-2 * ct.numst)));

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
double DFTsolver<OrbitalsType>::evaluateEnergy(
    const OrbitalsType& orbitals, const bool print_flag)
{
    // save energy recent history
    eks_history_[1] = eks_history_[0];

    // Get the new total energy
    const double ts = 0.5 * proj_matrices_->computeEntropy(); // in [Ha]
    eks_history_[0] = energy_->evaluateTotal(
        ts, proj_matrices_, ions_, orbitals, print_flag, os_);

    sum_eig_[1] = sum_eig_[0];
    sum_eig_[0] = 2. * proj_matrices_->getEigSum(); // 2.*sum in [Ry]

    return eks_history_[0];
}

template <class OrbitalsType>
int DFTsolver<OrbitalsType>::solve(OrbitalsType& orbitals,
    OrbitalsType& work_orbitals, Ions& ions, const short max_steps,
    const short iprint, double& last_eks)
{
    solve_tm_.start();

    Control& ct(*(Control::instance()));

    eks_history_[0] = 1.e9;
    eks_history_[1] = 1.e9;

    sum_eig_[0] = 1.e9;
    sum_eig_[1] = 1.e9;

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
        os_ << "### DFTsolver ###" << std::endl;
        if (ct.verbose > 1 && ct.precond_factor_computed)
            os_ << "Preconditioning factor: " << ct.precond_factor << std::endl;
    }

    orbitals_stepper_->setup(orbitals);

    if (ct.resetVH())
    {
        if (onpe0) os_ << "DFTsolver: reset Hartree potential" << std::endl;
        electrostat_->resetSolution();
    }

    // main electronic structure outer loop
    for (short step = 0; step <= max_steps; step++)
    {
        bool orthof = false;
        if (ct.orthof)
            orthof = ((((step + 1) % ct.orthof) == ct.max_changes_pot)
                      || ct.orthof == 1 || step == max_steps - 1);

        // turn on PB solver if necessary
        diel_control_.activate(it_scf_, deig2_);

        proj_matrices_->resetDotProductMatrices();

        // Generate new density
        rho_->update(orbitals);

        // Update potential
        if (testUpdatePot()) mgmol_strategy_->update_pot(ions);

        mgmol_strategy_->updateHmatrix(orbitals, ions);

        // theta = invB * Hij
        // (to be used for energy and gradient computation)
        proj_matrices_->updateThetaAndHB();

        if (step == max_steps) break;

        // Output the eigenvalues and occupations
        bool flag = false;
        if (iprint) flag = (!(step % iprint) && step < max_steps - 1);

        if (flag && !ct.short_sighted) mgmol_strategy_->printEigAndOcc();

        last_eks = evaluateEnergy(orbitals, flag);

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
                    << " DFTsolver: convergence achieved for delta E..."
                    << std::endl;
            break;
        }

        bool print_res = checkPrintResidual(step);

        const bool ortho
            = (ct.getOrthoType() == OrthoType::Eigenfunctions || orthof);

        if (!ortho || !ct.fullyOccupied())
        {
            // strip dm from the overlap contribution
            // dm <- Ls**T * dm * Ls
            dm_strategy_->stripDM();
        }

        // one step wave functions update
        // S and S^-1 should be up to date after that call
        const double restol = ct.checkResidual() ? ct.conv_tol : -1.;
        retval = orbitals_stepper_->updateWF(orbitals, ions, ct.precond_factor,
            orthof, work_orbitals, accelerate_, print_res, restol);

        if (ortho)
        {
            if (ct.isLocMode())
            {
                orbitals.normalize();
            }
            else
            {
                bool updateDM = false;
                bool updatedS = false;
                if (!ct.fullyOccupied())
                {
                    orbitals.computeGramAndInvS();
                    dm_strategy_->dressDM();
                    updateDM = true;
                    updatedS = true;
                }
                orbitals.orthonormalizeLoewdin(updatedS, nullptr, updateDM);

                orbitals_stepper_->restartMixing();
            }
        }
        else
        {
            orbitals.normalize();

            // if orthonorm() not called, recompute overlap
            orbitals.computeGramAndInvS();

            // rebuild dm with new overlap matrix
            dm_strategy_->dressDM();
        }

        if (retval == 0)
        {
            if (onpe0)
                os_ << std::endl
                    << std::endl
                    << " DFTsolver: convergence achieved for residual..."
                    << std::endl;
            break;
        }

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

        // updated Hij needed to compute new DM
        if (dm_strategy_->needH())
            mgmol_strategy_->updateHmatrix(orbitals, ions);

        // compute new density matrix
        dm_strategy_->update(orbitals);

        incInnerIt();

    } // end iterations

    if (iprint && !ct.short_sighted) mgmol_strategy_->printEigAndOcc();

#ifdef HAVE_ARPACK
    if (ct.precond_factor_computed)
    {
        const double alpha = 1. / getLAeigen(0., 500, ions);
        if (onpe0)
        {
            os_ << "Quench electrons" << std::endl;
            os_ << "A Posteriori Preconditioning factor: " << alpha
                << std::endl;
        }
    }
#endif
    solve_tm_.stop();

    return retval;
}

template <class OrbitalsType>
void DFTsolver<OrbitalsType>::printTimers(std::ostream& os)
{
    solve_tm_.print(os);
}

template class DFTsolver<LocGridOrbitals<ORBDTYPE>>;
template class DFTsolver<ExtendedGridOrbitals<ORBDTYPE>>;
