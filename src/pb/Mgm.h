// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PB_MGM_H
#define PB_MGM_H

#include "Vcycle.h"
#include "tools.h"
#include <iomanip>

namespace pb
{

template <class T1, class T2, typename T3>
bool Mgm(T1& A, T2& vh, const GridFunc<T3>& rho, const short cogr,
    const short max_sweeps, const double tol, const short nu1, const short nu2,
    const bool gather_coarse_level, double& final_residual,
    double& final_relative_residual, double& residual_reduction,
    short& nb_sweeps)
{
    // divide vh by inv(sqrt(epsilon)) (point by point)
    A.inv_transform(vh);

    const short bcx = vh.bc(0);
    const short bcy = vh.bc(1);
    const short bcz = vh.bc(2);

    const Grid& finegrid = vh.grid();

    // Compute r.h.s. from rho
    GridFunc<T3> res(rho);
    A.transform(res);
    GridFunc<T3> rhs(finegrid, bcx, bcy, bcz);
    A.rhs(res, rhs);

    // Hartree units
    // work GridFunc<T3>
    GridFunc<T3> lhs(finegrid, bcx, bcy, bcz);

    short bcwork[3] = { bcx, bcy, bcz };
    for (short d = 0; d < 3; d++)
        if (bcwork[d] == 2) bcwork[d] = 0;
    T2 work1(finegrid, bcwork[0], bcwork[1], bcwork[2]);

    const double inv_rhs_norm = 1. / norm(rhs);
    double init_residual_norm = 1.;
    bool converged            = false;
    nb_sweeps                 = 0;
    for (short i = 0; i < max_sweeps; i++)
    {

        A.apply(vh, lhs);
        // res=rhs-lhs;
        res.diff(rhs, lhs);

        const double res_norm = norm(res);
        if (i == 0) init_residual_norm = res_norm;

        // check reletive residual norm for convergence
        if (res_norm * inv_rhs_norm < tol)
        {
            final_residual          = res_norm;
            final_relative_residual = res_norm * inv_rhs_norm;
            converged               = true;
            break;
        }
#ifdef DEBUG
        GridFunc<T3> diff(finegrid, bcx, bcy, bcz);
        diff.diff(rhs, lhs);
        double resn = norm(diff);
        if (finegrid.mype_env().mytask() == 0)
        {
            cout << setprecision(3) << "Mgm: residual = " << resn << endl;
        }
#endif

        work1 = 0.;
        Vcycle(A, work1, res, cogr, nu1, nu2, gather_coarse_level);
        nb_sweeps++;

        vh += work1;
    }

    if (!converged)
    {
        A.apply(vh, lhs);
        lhs -= rhs;
        final_residual          = norm(lhs);
        final_relative_residual = final_residual * inv_rhs_norm;
    }

#ifdef DEBUG
    lhs.print_radial("residual.dat");
    if (finegrid.mype_env().mytask() == 0)
    {
        cout << setprecision(3) << "Mgm: final residual = " << final_residual
             << endl;
    }
#endif

    A.transform(vh);

    residual_reduction = final_residual / init_residual_norm;

    return converged;
}

} // namespace pb

#endif
