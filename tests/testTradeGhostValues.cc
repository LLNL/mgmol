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

#include <array>

// function of periodicity nx, ny, nz
double cos3(const int i, const int j, const int k, const int nx, const int ny,
    const int nz)
{
    double x = 2. * M_PI * (double)i / (double)nx;
    double y = 2. * M_PI * (double)j / (double)ny;
    double z = 2. * M_PI * (double)k / (double)nz;

    return std::cos(x) * std::cos(y) * std::cos(z);
}

TEST_CASE("Trade ghost values", "[trade]")
{
    const double origin[3]  = { 0., 0., 0. };
    const double ll         = 1.;
    const double lattice[3] = { ll, ll, ll };

    const int nfunc = 10;

    MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    // prepare 3 mesh sizes to test with
    std::vector<std::array<unsigned, 3>> meshes;
    {
        std::array<unsigned, 3> ngpts1 { { 32, 24, 20 } };
        meshes.push_back(ngpts1);

        std::array<unsigned, 3> ngpts2 { { 20, 32, 24 } };
        meshes.push_back(ngpts2);

        std::array<unsigned, 3> ngpts3 { { 24, 20, 32 } };
        meshes.push_back(ngpts3);
    }

    // test for several different number of mesh points
    // to generate different communications
    // patterns
    for (short nghosts = 1; nghosts < 3; nghosts++)
        for (auto& mesh : meshes)
        {
            // sync output
            mmpi.barrier();

            if (mmpi.PE0())
            {
                std::cout << "===================================" << std::endl;
                std::cout << "Number of ghosts points: " << nghosts
                          << std::endl;
            }

            const unsigned ngpts[3] = { mesh[0], mesh[1], mesh[2] };
            if (mmpi.PE0())
                std::cout << "Mesh " << ngpts[0] << " x " << ngpts[1] << " x "
                          << ngpts[2] << std::endl;

            pb::PEenv mype_env(MPI_COMM_WORLD, ngpts[0], ngpts[1], ngpts[2]);
            pb::Grid grid(origin, lattice, ngpts, mype_env, nghosts, 0);

            if (mmpi.PE0())
                std::cout << "MPI tasks grid is " << mype_env.n_mpi_task(0)
                          << " x " << mype_env.n_mpi_task(1) << " x "
                          << mype_env.n_mpi_task(2) << std::endl;

            if (mmpi.PE0()) std::cout << "GridFunc..." << std::endl;

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
                        inner_data[iiy + iz]
                            = cos3(ix, iy, iz, nx, ny, nz) + uvalue;
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
                        CHECK(uu[iiy + iz]
                              == Approx(cos3(ix - nghosts, iy - nghosts,
                                            iz - nghosts, nx, ny, nz)
                                        + uvalue)
                                     .epsilon(1.e-8));
                    }
                }
            }

            mmpi.barrier();
            if (mmpi.PE0()) std::cout << "GridFuncVector..." << std::endl;

            // now test with GridFuncVector
            std::vector<std::vector<int>> gids;
            gids.resize(1);
            for (int i = 0; i < nfunc; i++)
                gids[0].push_back((i + 3) % nfunc);

            pb::GridFuncVector<double> gfv(grid, 1, 1, 1, gids);
            for (int i = 0; i < nfunc; i++)
            {
                std::vector<double> scaled_data(inner_data);
                LinearAlgebraUtils<MemorySpace::Host>::MPscal(
                    nx * ny * nz, (double)(i + 1), scaled_data.data());
                gfv.assign(i, scaled_data.data(), 'd');
            }
            gfv.trade_boundaries();

            if (mype_env.onpe0())
                for (int i = 0; i < nfunc; i++)
                {
                    const pb::GridFunc<double>& gfi(gfv.getGridFunc(i));
                    double* uu = gfi.uu();
                    for (int ix = initx; ix < endx; ix++)
                    {
                        int iix = ix * grid.inc(0);

                        for (int iy = inity; iy < endy; iy++)
                        {
                            int iiy = iy * grid.inc(1) + iix;

                            for (int iz = initz; iz < endz; iz++)
                            {
                                double ref_val
                                    = (i + 1)
                                      * (uvalue
                                          + cos3(ix - nghosts, iy - nghosts,
                                              iz - nghosts, nx, ny, nz));

                                CHECK(uu[iiy + iz]
                                      == Approx(ref_val).epsilon(1.e-8));
                            }
                        }
                    }
                }
        }
}
