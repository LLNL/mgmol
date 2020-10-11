// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "GridFuncVector.h"
#include "MGmol_MPI.h"
#include "PEenv.h"
#include "mputils.h"

#include "catch.hpp"

TEST_CASE("Trade ghost values", "[trade]")
{
    const double origin[3]  = { 0., 0., 0. };
    const double ll         = 1.;
    const double lattice[3] = { ll, ll, ll };
    const short nghosts     = 1;

    const int nfunc = 10;

    MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

    // test for several different number of mesh points
    // in z-direction to generate different communications
    // patterns
    for (unsigned nzg = 20; nzg <= 80; nzg *= 2)
    {
        const unsigned ngpts[3] = { 32, 24, nzg };
        std::cout << "Mesh " << ngpts[0] << " x " << ngpts[1] << " x "
                  << ngpts[2] << std::endl;
        pb::PEenv mype_env(MPI_COMM_WORLD, ngpts[0], ngpts[1], ngpts[2]);
        pb::Grid grid(origin, lattice, ngpts, mype_env, nghosts, 0);

        std::cout << "MPI tasks grid is " << mype_env.n_mpi_task(0) << " x "
                  << mype_env.n_mpi_task(1) << " x " << mype_env.n_mpi_task(2)
                  << std::endl;

        // periodic GridFunc
        pb::GridFunc<double> gf(grid, 1, 1, 1);

        int nx = static_cast<int>(grid.dim(0));
        int ny = static_cast<int>(grid.dim(1));
        int nz = static_cast<int>(grid.dim(2));

        // set interior values to uniform value
        std::vector<double> inner_data(nx * ny * nz);
        const double uvalue = 1.111;
        for (int ix = 0; ix < nx; ix++)
        {
            int iix = ix * ny * nz;

            for (int iy = 0; iy < ny; iy++)
            {
                int iiy = iy * nz + iix;

                for (int iz = 0; iz < nz; iz++)
                {
                    inner_data[iiy + iz] = uvalue;
                }
            }
        }

        // initialize gf with inner_data
        gf.assign(inner_data.data(), 'd');
        gf.set_updated_boundaries(false);

        double norm_before = gf.norm2();

        // fill ghost values with data from neighboring subdomain
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

        double* uu = gf.uu();
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

        // now test with GridFuncVector
        std::vector<std::vector<int>> gids;
        gids.resize(1);
        for (int i = 0; i < nfunc; i++)
            // gids[0].push_back((i + 3) % nfunc);
            gids[0].push_back(i);

#ifndef HAVE_OPENMP_OFFLOAD
        pb::GridFuncVector<double, MemorySpace::Host> gfv(grid, 1, 1, 1, gids);
        for (int i = 0; i < nfunc; i++)
        {
            std::vector<double> scaled_data(inner_data);
            LinearAlgebraUtils<MemorySpace::Host>::MPscal(
                nx * ny * nz, (double)(i + 1), scaled_data.data());
            gfv.assign(i, scaled_data.data(), 'd');
        }
        gfv.trade_boundaries();
#else
        pb::GridFuncVector<double, MemorySpace::Device> gfv(
            grid, 1, 1, 1, gids);
        for (int i = 0; i < nfunc; i++)
        {
            std::vector<double> scaled_data(inner_data);
            LinearAlgebraUtils<MemorySpace::Host>::MPscal(
                nx * ny * nz, (double)(i + 1), scaled_data.data());
            gfv.assign(i, scaled_data.data(), 'd');
        }
        gfv.copyHtoD(nfunc * grid.sizeg());
        gfv.trade_boundaries();
        gfv.copyDtoH(nfunc * grid.sizeg());
#endif

        if (mype_env.onpe0())
            for (int i = 0; i < nfunc; i++)
            {
                const pb::GridFunc<double>& gfi(gfv.getGridFunc(i));
                double* uu     = gfi.uu();
                double ref_val = (i + 1) * uvalue;
                for (int ix = initx; ix < endx; ix++)
                {
                    int iix = ix * grid.inc(0);

                    for (int iy = inity; iy < endy; iy++)
                    {
                        int iiy = iy * grid.inc(1) + iix;

                        for (int iz = initz; iz < endz; iz++)
                        {
                            CHECK(
                                uu[iiy + iz] == Approx(ref_val).epsilon(1.e-6));
                        }
                    }
                }
            }
    }
}
