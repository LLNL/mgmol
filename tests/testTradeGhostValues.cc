// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "GridFunc.h"
#include "MGmol_MPI.h"
#include "PEenv.h"

#include "catch.hpp"

TEST_CASE("Trade ghost values", "[trade]")
{
    const double origin[3]  = { 0., 0., 0. };
    const double ll         = 1.;
    const double lattice[3] = { ll, ll, ll };
    const unsigned ngpts[3] = { 32, 24, 20 };
    const short nghosts     = 1;

    MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

    pb::PEenv mype_env(MPI_COMM_WORLD, ngpts[0], ngpts[1], ngpts[2]);
    pb::Grid grid(origin, lattice, ngpts, mype_env, nghosts, 0);
    // periodic GridFunc
    pb::GridFunc<double> gf(grid, 1, 1, 1);

    // set interior values to uniform value
    double* uu          = gf.uu();
    const double uvalue = 1.111;
    for (int ix = nghosts; ix < static_cast<int>(nghosts + grid.dim(0)); ix++)
    {
        int iix = ix * grid.inc(0);

        for (int iy = nghosts; iy < static_cast<int>(nghosts + grid.dim(1));
             iy++)
        {
            int iiy = iy * grid.inc(1) + iix;

            for (int iz = nghosts; iz < static_cast<int>(nghosts + grid.dim(2));
                 iz++)
            {
                uu[iiy + iz] = uvalue;
            }
        }
    }

    gf.set_updated_boundaries(false);

    double norm_before = gf.norm2();

    gf.trade_boundaries();

    // check norm has not changed
    double norm_after = gf.norm2();
    CHECK(norm_after == Approx(norm_before).epsilon(1.e-8));

    // check ghost values
    const int initx = 0;
    const int inity = 0;
    const int initz = 0;

    const int endx = 2 * nghosts + grid.dim(0);
    const int endy = 2 * nghosts + grid.dim(1);
    const int endz = 2 * nghosts + grid.dim(2);

    for (int ix = initx; ix < endx; ix++)
    {
        int iix = ix * grid.inc(0);

        for (int iy = inity; iy < endy; iy++)
        {
            int iiy = iy * grid.inc(1) + iix;

            for (int iz = initz; iz < endz; iz++)
            {
                CHECK(uu[iiy + iz] == Approx(uvalue).epsilon(1.e-6));
            }
        }
    }
}
