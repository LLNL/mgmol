// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ABPG.h"
#include "Control.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "OrthoAndersonMix.h"
#include "Potentials.h"
#include "ProjectedMatricesInterface.h"

#include <iostream>

template <class T>
void ABPG<T>::setup(T& orbitals)
{
    Control& ct = *(Control::instance());

    if (ct.wf_dyn == 1) // use Anderson extrapolation
    {
        if (ct.getOrthoType() == OrthoType::Orthonormal)
            wf_mix_
                = new OrthoAndersonMix<T>(ct.wf_m, ct.betaAnderson, orbitals);
        else
            wf_mix_ = new AndersonMix<T>(ct.wf_m, ct.betaAnderson, orbitals);
    }
}

//
// Performs a single wave functions update step.
//
// orthof=true: wants orthonormalized updated wave functions
template <class T>
int ABPG<T>::updateWF(T& orbitals, Ions& ions, const double precond_factor,
    const bool /*orthof*/, T& work_orbitals, const bool accelerate,
    const bool print_res, const double atol)
{
    abpg_nl_update_tm_.start();

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2) os_ << "ABPG::update()..." << std::endl;

    // temporary Orbitals to hold residual
    T res("Residual", orbitals, false);

    const bool check_res = (atol > 0.);
    double normRes = mgmol_strategy_->computeResidual(orbitals, work_orbitals,
        ions, res, (print_res || check_res), check_res);
    if (normRes < atol && check_res)
    {
        abpg_nl_update_tm_.stop();
        return 0;
    }

    // update wavefunctions
    update_states(orbitals, res, work_orbitals, precond_factor, accelerate);

    abpg_nl_update_tm_.stop();

    return 1;
}

//////////////////////////////////////////////////////////////////////////////
// Update orbitals using MG preconditioning and Anderson extrapolation
template <class T>
void ABPG<T>::update_states(T& orbitals, T& res, T& work_orbitals,
    const double precond_factor, const bool accelerate)
{
    assert(orbitals.getIterativeIndex() >= 0);

    update_states_tm_.start();

    Control& ct = *(Control::instance());

    if (ct.withPreconditioner())
    {
        // PRECONDITIONING
        // compute the preconditioned steepest descent direction
        // -> res
        mgmol_strategy_->precond_mg(res);
    }

    // non-energy compatible spread penalties are added after
    // preconditioning is applied
    if (ct.isSpreadFunctionalActive() && !ct.isSpreadFunctionalEnergy())
        mgmol_strategy_->addResidualSpreadPenalty(orbitals, res);

    // apply AOMM projector to have a gradient orthogonal to kernel functions
    if (ct.use_kernel_functions)
    {
        mgmol_strategy_->projectOutKernel(res);

        // just reapply mask for now...
        res.applyMask();
    }

    if (ct.project_out_psd)
    {
        if (onpe0) os_ << "Project out preconditioned gradient" << std::endl;
        orbitals.projectOut(res);
    }

    const double alpha = 0.5 * precond_factor;

    // Update wavefuntion
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "Update states" << std::endl;
#else
    if (onpe0 && ct.verbose > 2) os_ << "Update states" << std::endl;
#endif
    if (accelerate)
    {
        res.scal(alpha);
        // Extrapolation scheme
        assert(wf_mix_ != nullptr);
        std::ostream os(nullptr);
        if (onpe0) os.rdbuf(std::cout.rdbuf());
        wf_mix_->update(res, work_orbitals, os, (ct.verbose > 0));
    }
    else
    {
        // Preconditioned Power Method
        orbitals.axpy((ORBDTYPE)alpha, res);

        if (ct.getOrthoType() == OrthoType::Orthonormal)
            orbitals.orthonormalizeLoewdin(false);
    }
    orbitals.incrementIterativeIndex();

    update_states_tm_.stop();

    assert(orbitals.getIterativeIndex() >= 0);
}

template <class T>
void ABPG<T>::printTimers(std::ostream& os)
{
    abpg_tm_.print(os);
    abpg_nl_update_tm_.print(os);
    comp_res_tm_.print(os);
    update_states_tm_.print(os);
}

template class ABPG<LocGridOrbitals<ORBDTYPE>>;
template class ABPG<ExtendedGridOrbitals<ORBDTYPE>>;
