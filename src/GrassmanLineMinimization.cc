// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "MGmol.h"

#include "Control.h"
#include "GrassmanLineMinimization.h"
#include "Hamiltonian.h"
#include "LocGridOrbitals.h"
#include "Potentials.h"

template <class T>
Timer GrassmanLineMinimization<T>::line_min_tm_("Grassman_line_min");
template <class T>
Timer GrassmanLineMinimization<T>::nl_update_tm_("Grassman_nl_update");
template <class T>
Timer GrassmanLineMinimization<T>::comp_res_tm_("Grassman_comp_res");
template <class T>
Timer GrassmanLineMinimization<T>::update_states_tm_("Grassman_update_states");

//
// Performs a single self consistent step.
//
// orthof=true: wants orthonormalized updated wave functions
template <class T>
int GrassmanLineMinimization<T>::updateWF(T& orbitals, Ions& ions,
    const double precond_factor, const bool orthof, T& work_orbitals,
    const bool accelerate, const bool print_res, const double atol)
{
    nl_update_tm_.start();

    (void)orthof; // not used

    static bool first_time = true;

    if (first_time)
    {
        first_time = false;
        conjugate_ = false;

        new_grad_   = new T("NewG", orbitals, false);
        new_pcgrad_ = new T("NewP", *new_grad_);
    }
    else
    {
        conjugate_ = true;
    }

    accelerate_ = accelerate;

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
        os_ << "GrassmanLineMinimization<T>::update() ..." << std::endl;

    // Update wavefunctions
    const bool check_res = (atol > 0.);
    double normRes = mgmol_strategy_->computeResidual(orbitals, work_orbitals,
        ions, *new_grad_, (print_res || check_res), check_res);
    if (normRes < atol && check_res)
    {
        nl_update_tm_.stop();
        return 0;
    }

    // update states
    update_states(orbitals, precond_factor);

    nl_update_tm_.stop();

    return 1;
}

//////////////////////////////////////////////////////////////////////////////

template <class T>
void GrassmanLineMinimization<T>::update_states(
    T& orbitals, const double precond_factor)
{
    assert(orbitals.getIterativeIndex() >= 0);

    update_states_tm_.start();

    Control& ct = *(Control::instance());

    // compute preconditioned residual
    new_pcgrad_->assign(*new_grad_);
    if (ct.withPreconditioner())
    {
        // PRECONDITIONING
        mgmol_strategy_->precond_mg(*new_pcgrad_);
    }
    const double alpha = 0.5 * precond_factor;

    // Update wavefuntion
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "Update states" << std::endl;
#else
    if (onpe0 && ct.verbose > 2) os_ << "Update states" << std::endl;
#endif

    // do conjugation
    conjugate();
    // compute step size and update functions
    if (accelerate_)
    {
        // Grassman line minimization method
        double lambda = computeStepSize(orbitals);
        //        orbitals.projectOut(*sdir_);
        orbitals.axpy((ORBDTYPE)lambda, *sdir_);
        // recompute overlap and inverse for new wavefunctions
        orbitals.computeGramAndInvS();
        if (onpe0 && ct.verbose > 1)
            std::cout << "Grassman CG: lambda = " << lambda << std::endl;
        // now do parallel transport update of history data
        parallelTransportUpdate(lambda, orbitals);
    }
    else
    {
        // Preconditioned Power Method
        //        orbitals.projectOut(*sdir_);
        orbitals.axpy((ORBDTYPE)alpha, *sdir_);
        //        if(onpe0)cout<<"alpha = "<<alpha<<endl;
        // recompute overlap and inverse for new wavefunctions
        orbitals.computeGramAndInvS();
        // now do parallel transport update of history data
        parallelTransportUpdate(alpha, orbitals);
    }
    orbitals.incrementIterativeIndex();

    update_states_tm_.stop();

    assert(orbitals.getIterativeIndex() >= 0);
}

template <class T>
void GrassmanLineMinimization<T>::printTimers(std::ostream& os)
{
    line_min_tm_.print(os);
    nl_update_tm_.print(os);
    comp_res_tm_.print(os);
    update_states_tm_.print(os);
}

template class GrassmanLineMinimization<LocGridOrbitals<ORBDTYPE>>;
template class GrassmanLineMinimization<ExtendedGridOrbitals<ORBDTYPE>>;
