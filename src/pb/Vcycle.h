// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PB_VCYCLE_H
#define PB_VCYCLE_H

#define VCYCLE_DEBUG 0
#define USE_LOWER_ORDER 1 // decrease operator order when coarsening

#include "GridFunc.h"
#include <iostream>

extern Timer vcycle_repwrk_tm;
extern Timer vcycle_repinit_tm;
extern Timer vcycle_repend_tm;
extern Timer vcycle_repvcycle_tm;

namespace pb
{

// assumes x=0 in input
template <class T1, class T2, typename T3>
int Vcycle(T1& A, T2& x, const GridFunc<T3>& rhs, const short cogr,
    const short nu1, const short nu2, const bool gather_coarse_level)
{
    assert(x.grid().dim(0) < 10000);
    assert(x.grid().dim(1) < 10000);
    assert(x.grid().dim(2) < 10000);

    const Grid& level_grid = A.grid();
    const short nghosts    = x.ghost_pt();

#if VCYCLE_DEBUG
    if (level_grid.mype_env().onpe0())
        std::cout << " Vcycle at level " << level_grid.level()
                  << ", dim0=" << level_grid.dim(0) << ", nghosts=" << nghosts
                  << ", cogr=" << cogr << std::flush << std::endl;
    level_grid.mype_env().barrier();

    double norm_tmp = norm(rhs);
    if (level_grid.mype_env().onpe0())
        std::cout << "Vcycle: rhs=" << norm_tmp << std::endl;
    norm_tmp = norm(x);
    if (level_grid.mype_env().onpe0())
        std::cout << "Vcycle: x=" << norm_tmp << std::endl;
#endif

    const bool flag_coarsen
        = ((!(level_grid.dim(0)
               & 1)) // cannot coarsen if mesh not divisible by 2
            && (!(level_grid.dim(1) & 1)) && (!(level_grid.dim(2) & 1))
            && (static_cast<int>(level_grid.dim(0)) >= 2 * nghosts)
            && (static_cast<int>(level_grid.dim(1)) >= 2 * nghosts)
            && (static_cast<int>(level_grid.dim(2)) >= 2 * nghosts));
#if VCYCLE_DEBUG
    if (level_grid.mype_env().onpe0())
        std::cout << "Vcycle: flag_coarsen=" << flag_coarsen << std::endl;
#endif

    const bool flag_gather = ((!flag_coarsen || (level_grid.level() <= (-cogr)))
                              && gather_coarse_level);

    // try gather data on single processor to reach coarser levels
    if (flag_gather && level_grid.mype_env().n_mpi_tasks() > 1)
    {
        vcycle_repwrk_tm.start();
#if VCYCLE_DEBUG
        if (level_grid.mype_env().onpe0()) cout << "Vcycle: gather" << endl;
#endif

        vcycle_repinit_tm.start();
        const PEenv replicated_peenv(level_grid.mype_env().comm_global());
        Grid replicated_grid(level_grid.replicated_grid(replicated_peenv));

#if 1 // gather and solve on all PEs
        T3* replicated_func = new T3[level_grid.gsize()];

        // gather rhs
        rhs.init_vect(replicated_func, 'g');
        GridFunc<T3> replicated_rhs(replicated_grid, x.bc(0), x.bc(1), x.bc(2));
        replicated_rhs.assign(replicated_func);

        GridFunc<T3> replicated_x(replicated_grid, x.bc(0), x.bc(1), x.bc(2));

        T1 replicated_A = A.replicatedOp(replicated_grid);

        vcycle_repinit_tm.stop();
        vcycle_repvcycle_tm.start();
        // solve
        int ret = Vcycle(replicated_A, replicated_x, replicated_rhs,
            cogr - level_grid.level(), nu1, nu2, false);
        vcycle_repvcycle_tm.stop();

        vcycle_repend_tm.start();
        // scatter
        x.assign(replicated_x, 'g');

        delete[] replicated_func;

#else // gather and solve on PE 0

        const bool onpe0    = level_grid.mype_env().onpe0();
        T3* replicated_func = NULL;
        if (onpe0)
        {
            replicated_func = new T3[level_grid.gsize()];
        }
        rhs.gather(replicated_func);
        GridFunc<T3> replicated_x(replicated_grid, x.bc(0), x.bc(1), x.bc(2));

        int ret = -1;
        if (onpe0)
        {
            GridFunc<T3> replicated_rhs(
                replicated_func, replicated_grid, x.bc(0), x.bc(1), x.bc(2), 0);

            T1 replicated_A = A.replicatedOp(replicated_grid);

            vcycle_repinit_tm.stop();
            vcycle_repvcycle_tm.start();
            // solve
            ret = Vcycle(replicated_A, replicated_x, replicated_rhs,
                cogr - level_grid.level(), nu1, nu2, false);
            vcycle_repvcycle_tm.stop();

            vcycle_repend_tm.start();
            delete[] replicated_func;
        }

        // scatter
        x.scatterFrom(replicated_x);
#endif

        vcycle_repend_tm.stop();
        vcycle_repwrk_tm.stop();

        return ret;
    }

    GridFunc<T3> res(x.grid(), x.bc(0), x.bc(1), x.bc(2));

    // pre-smoothing
    for (short i = 0; i < nu1; i++)
    {

        A.jacobi(x, rhs, res);
#if VCYCLE_DEBUG
        norm_tmp = norm(res);
        if (level_grid.mype_env().onpe0())
            std::cout << "pre-smoothing level " << level_grid.level()
                      << ", residual=" << norm_tmp << std::endl;
#endif
    }

    if (level_grid.level() > (-cogr))
    {
        if (flag_coarsen)
        {

            A.apply(x, res);
            res -= rhs; // get residual
#if VCYCLE_DEBUG
            norm_tmp = norm(res);
            if (level_grid.mype_env().onpe0())
                cout << "level " << level_grid.level()
                     << ", residual=" << norm_tmp << endl;
#endif

            // restrictions
            T1 B = A.coarseOp(level_grid);
            // GridFunc<T3>  rcoarse(B.grid(),x.bc(0),x.bc(1),x.bc(2));
#if USE_LOWER_ORDER
            B.setLowerOrderGrid();
            int ng = std::max(level_grid.ghost_pt() - 1, 1);
            const Grid coarse_grid(B.getLowerOrderGrid(), ng);
#else
            const Grid coarse_grid(B.grid());
#endif
            assert(coarse_grid.dim(0) < 10000);
            assert(coarse_grid.dim(1) < 10000);
            assert(coarse_grid.dim(2) < 10000);

            GridFunc<T3> rcoarse(coarse_grid, x.bc(0), x.bc(1), x.bc(2));
            res.restrict3D(rcoarse);

            short bc[3] = { x.bc(0), x.bc(1), x.bc(2) };
            for (short d = 0; d < 3; d++)
                if (bc[d] == 2) bc[d] = 0;

            T2 ucoarse(coarse_grid, bc[0], bc[1], bc[2]);

#if USE_LOWER_ORDER
            Vcycle(B.getLowerOrderOp(), ucoarse, rcoarse, cogr, nu1, nu2,
                gather_coarse_level);
#else
            Vcycle(B, ucoarse, rcoarse, cogr, nu1, nu2, gather_coarse_level);
#endif

            res.extend3D(ucoarse);
            x -= res;
        }
        else
        {
#if VCYCLE_DEBUG
            level_grid.mype_env().barrier();
            if (level_grid.mype_env().onpe0())
                std::cout << " Coarsest level: mesh " << x.dim(0) << "x"
                          << x.dim(1) << "x" << x.dim(2) << " reached\n";
#endif
        }
    }

    // post-smoothing
    for (short i = 0; i < nu2; i++)
    {

        A.jacobi(x, rhs, res);
        // if(x.bc()>0)x.average0();
#if VCYCLE_DEBUG
        norm_tmp = norm(res);
        if (level_grid.mype_env().onpe0())
            std::cout << "post-smoothing level " << level_grid.level()
                      << ", residual=" << norm_tmp << std::endl;
#endif
    }
#if VCYCLE_DEBUG
    norm_tmp = norm(res);
    if (level_grid.mype_env().onpe0())
        if (level_grid.level() == 0)
            std::cout << "level " << level_grid.level()
                      << ", residual=" << norm_tmp << std::endl;
    norm_tmp = norm(x);
    if (level_grid.mype_env().onpe0())
        cout << "end of Vcycle: norm(x)=" << norm_tmp << endl;
    A.apply(x, res);
    res -= rhs;
    norm_tmp = norm(res);
    if (level_grid.mype_env().onpe0())
        std::cout << "end of Vcycle: norm(res)=" << norm_tmp << std::endl;
#endif

    return 0;
}

} // namespace pb

#endif
