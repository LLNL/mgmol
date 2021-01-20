#include "GridFuncVector.h"
#include "PEenv.h"

#include "catch.hpp"

#include <vector>

TEST_CASE("MG kernets", "[MGkernels]")
{
    const double origin[3]  = { 0., 0., 0. };
    const double ll         = 2.;
    const double lattice[3] = { ll, ll, ll };
    const unsigned ngpts[3] = { 32, 24, 20 };
    const short nghosts     = 2;

    const int nfunc = 5;

    pb::PEenv mype_env(MPI_COMM_WORLD, ngpts[0], ngpts[1], ngpts[2]);
    pb::Grid grid(origin, lattice, ngpts, mype_env, nghosts, 0);

    // periodic GridFuncVector
    std::vector<std::vector<int>> gids;
    gids.resize(1);
    for (int i = 0; i < nfunc; i++)
        gids[0].push_back(i);

    pb::GridFuncVector<double> gffine(grid, 1, 1, 1, gids);

    const double value = 1.11;

    // initialize gffine without ghost values
    double* ufine  = gffine.data();
    const int endx = nghosts + grid.dim(0);
    const int endy = nghosts + grid.dim(1);
    const int endz = nghosts + grid.dim(2);

    for (int ifunc = 0; ifunc < nfunc; ifunc++)
    {
        int offset = ifunc * grid.sizeg();
        for (int ix = nghosts; ix < endx; ix++)
        {
            int iix = ix * grid.inc(0) + offset;
            for (int iy = nghosts; iy < endy; iy++)
            {
                int iiy = iy * grid.inc(1) + iix;
                for (int iz = nghosts; iz < endz; iz++)
                {
                    ufine[iiy + iz] = ifunc * value;
                }
            }
        }
    }
    gffine.set_updated_boundaries(false);

    pb::Grid coarse_grid(grid.coarse_grid());
    pb::GridFuncVector<double> gfcoarse(coarse_grid, 1, 1, 1, gids);

    // coarsen data
    gffine.restrict3D(gfcoarse);

    // refine data just coarsened into new object
    pb::GridFuncVector<double> gffine2(grid, 1, 1, 1, gids);
    gffine2.extend3D(gfcoarse);

    // check extended values
    double* ufine2 = gffine2.data();
    for (int ifunc = 0; ifunc < nfunc; ifunc++)
    {
        int offset = ifunc * grid.sizeg();
        for (int ix = nghosts; ix < endx; ix++)
        {
            int iix = ix * grid.inc(0) + offset;
            for (int iy = nghosts; iy < endy; iy++)
            {
                int iiy = iy * grid.inc(1) + iix;
                for (int iz = nghosts; iz < endz; iz++)
                {
                    double expected_val = ifunc * value;
                    CHECK(
                        ufine2[iiy + iz] == Approx(expected_val).margin(1.e-8));
                }
            }
        }
    }
}
