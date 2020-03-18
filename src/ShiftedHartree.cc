// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include <iomanip>
#include <iostream>
using namespace std;

#include "ShiftedHartree.h"

template class ShiftedHartree<pb::ShiftedLaph4M<POTDTYPE>>;
// template class ShiftedHartree<pb::ShiftedLaph4M<POTDTYPE> >;

// Solve Poisson problem
// For ct.screening_const>0., solve the problem with a shift according to
// Ref. M. Manninen et al., PRB 12, no 10, p. 4012 (see Appendix)
// Known as Kerker mixing in PW

template <class T>
void ShiftedHartree<T>::solve(
    const pb::GridFunc<RHODTYPE>& rho, const pb::GridFunc<RHODTYPE>& rhoc)
{
    PoissonInterface::poisson_tm_.start();

    pb::GridFunc<RHODTYPE> work_rho(rho);

    // Keep in memory vh*rho before updating vh
    const double vel        = Poisson::grid_.vel();
    Poisson::Int_vhrho_old_ = vel * Poisson::vh_->gdot(rho);

    // Subtract compensating charges from rho
    work_rho -= rhoc;

    work_rho.axpy(shift_ / (4. * M_PI), *Poisson::vh_);

    /* Check for uniform precision before calling poisson_solver.
     * Downgrade or upgrade rhs (work_rho) to have precision of solution (vh_).
     * Note that this could be done at the beginning of this function, but
     * several operations involving rho might be done in lower precision
     * (depending on POTDTYPE), which could affect accuracy. For now, we delay
     * the switch until just before the solve call.
     */
    //    if(sizeof(POTDTYPE) != sizeof(RHODTYPE))
    //    {
    /* solve with POTDTYPE precision */
    pb::GridFunc<POTDTYPE> rhs(work_rho);
    poisson_solver_->solve(*Poisson::vh_, rhs);
    //    }
    //    else
    //    {
    //       poisson_solver_->solve(*Poisson::vh_, work_rho);
    //    }

    double residual_reduction = poisson_solver_->getResidualReduction();
    double final_residual     = poisson_solver_->getFinalResidual();

    if (onpe0)
        (*MPIdata::sout) << setprecision(2) << scientific
                         << " Hartree: residual reduction = "
                         << residual_reduction
                         << ", final residual = " << final_residual << endl;

    Poisson::Int_vhrho_  = vel * Poisson::vh_->gdot(rho);
    Poisson::Int_vhrhoc_ = vel * Poisson::vh_->gdot(rhoc);

    PoissonInterface::poisson_tm_.stop();
}
