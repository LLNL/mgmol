// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <iomanip>
#include <iostream>
using namespace std;

#include "Control.h"
#include "Hartree_CG.h"
#include "MultipoleExpansion.h"

#include "Laph2.h"
#include "Laph4.h"
#include "Laph4M.h"
#include "Laph4MP.h"
#include "Laph6.h"
#include "Laph8.h"

// Timer Poisson::poisson_tm_("Poisson::poisson");

template <class T>
void Hartree_CG<T>::solve(
    const pb::GridFunc<RHODTYPE>& rho, const pb::GridFunc<RHODTYPE>& rhoc)
{
    PoissonInterface::poisson_tm_.start();

    pb::GridFunc<RHODTYPE> work_rho(rho);
    Control& ct = *(Control::instance());

    // Keep in memory vh*rho before updating vh
    const double vel        = Poisson::grid_.vel();
    Poisson::Int_vhrho_old_ = vel * Poisson::vh_->gdot(rho);

    //(*MPIdata::sout)<<"Integral rho="<<rho.integral()<<endl;
    //(*MPIdata::sout)<<"Integral rhoc="<<rhoc.integral()<<endl;

    // Subtract compensating charges from rho
    work_rho -= rhoc;

    int dim_mpol = 0;
    for (int i = 0; i < 3; i++)
        if (Poisson::bc_[i] == 2) dim_mpol++;
    //(*MPIdata::sout)<<"dim_mpol="<<dim_mpol<<endl;

    pb::GridFunc<POTDTYPE> bc_func(
        Poisson::grid_, Poisson::bc_[0], Poisson::bc_[1], Poisson::bc_[2]);
    if (dim_mpol > 0)
    {
        const Vector3D origin_cell(Poisson::grid_.origin(0),
            Poisson::grid_.origin(1), Poisson::grid_.origin(2));
        const Vector3D cell(
            Poisson::grid_.ll(0), Poisson::grid_.ll(1), Poisson::grid_.ll(2));

        MultipoleExpansion mp(Poisson::grid_, Poisson::bc_, origin_cell, cell);
        mp.setOrder(ct.multipole_order);
        mp.setup(work_rho);

        if (dim_mpol == 2) mp.expand2d(bc_func);
        if (dim_mpol == 3) mp.expand(bc_func);
        if (dim_mpol > 0)
        {
            Poisson::vh_->set_bc_func(&bc_func);
            // bc_func.print_radial("bcfunc");
        }
    }

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
    rhs *= (4. * M_PI);
    poisson_solver_->solve(*Poisson::vh_, rhs);
    //    }
    //    else
    //    {
    //       poisson_solver_->solve(*Poisson::vh_, work_rho);
    //    }

    double residual_reduction = poisson_solver_->getResidualReduction();
    double final_residual     = poisson_solver_->getFinalResidual();
    const bool large_residual
        = (residual_reduction > 1.e-3 || final_residual > 1.e-3);

    if (onpe0 && (large_residual || ct.verbose > 1))
        (*MPIdata::sout) << setprecision(2) << scientific
                         << "Hartree_CG: residual reduction = "
                         << residual_reduction
                         << ", final residual = " << final_residual << endl;

    Poisson::Int_vhrho_  = vel * Poisson::vh_->gdot(rho);
    Poisson::Int_vhrhoc_ = vel * Poisson::vh_->gdot(rhoc);

    PoissonInterface::poisson_tm_.stop();

    assert(residual_reduction == residual_reduction);
}

template class Hartree_CG<pb::Laph2<POTDTYPE>>;
// template class Hartree_CG<pb::Laph2<float> >;
template class Hartree_CG<pb::Laph4<POTDTYPE>>;
// template class Hartree_CG<pb::Laph4<float> >;
template class Hartree_CG<pb::Laph6<POTDTYPE>>;
// template class Hartree_CG<pb::Laph6<float> >;
template class Hartree_CG<pb::Laph8<POTDTYPE>>;
// template class Hartree_CG<pb::Laph8<float> >;
template class Hartree_CG<pb::Laph4M<POTDTYPE>>;
// template class Hartree_CG<pb::Laph4M<float> >;
template class Hartree_CG<pb::Laph4MP<POTDTYPE>>;
// template class Hartree_CG<pb::Laph4MP<float> >;
