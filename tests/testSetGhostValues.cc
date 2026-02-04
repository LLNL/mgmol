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

#include "catch.hpp"

#include <array>

TEST_CASE("Set ghost values", "[set ghosts")
{
    const double origin[3]  = { 0., 0., 0. };
    const double ll         = 1.;
    const double lattice[3] = { ll, ll, ll };

    MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    // prepare 3 mesh sizes to test with
    std::vector<std::array<unsigned, 3>> meshes;
    {
        std::array<unsigned, 3> ngpts1 { { 16, 24, 20 } };
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

            // non-periodic GridFunc
            pb::GridFunc<double> gf(grid, 0, 0, 0);

            int nx = static_cast<int>(grid.dim(0));
            int ny = static_cast<int>(grid.dim(1));
            int nz = static_cast<int>(grid.dim(2));

            if (mmpi.PE0())
                std::cout << "Local Mesh " << nx << " x " << ny << " x " << nz
                          << std::endl;

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

            // fill ghost values
            const bool direction[3] = { true, true, true };
            gf.defaultTrade_boundaries();
            gf.setBoundaryValues(uvalue, direction);

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
                        CHECK(uu[iiy + iz] == Approx(uvalue).epsilon(1.e-8));
                    }
                }
            }
        }
}
