#include "GridFunc.h"
#include "Laph4.h"
#include "PEenv.h"

#include "catch.hpp"

TEST_CASE("Laplacian 4th order", "[laph4]")
{
    const double origin[3]  = { 0., 0., 0. };
    const double ll         = 2.;
    const double lattice[3] = { ll, ll, ll };
    const unsigned ngpts[3] = { 32, 24, 20 };
    const short nghosts     = 2;

    const double h[3] = { ll / (double(ngpts[0])), ll / (double(ngpts[1])),
        ll / (double(ngpts[2])) };
    pb::PEenv mype_env(MPI_COMM_WORLD, ngpts[0], ngpts[1], ngpts[2]);
    pb::Grid grid(origin, lattice, ngpts, mype_env, nghosts, 0);

    pb::Laph4<double> lap(grid);

    // periodic GridFunc
    pb::GridFunc<double> gf1(grid, 1, 1, 1);
    pb::GridFunc<double> gf2(grid, 1, 1, 1);

    // initialize gf1
    double* u1     = gf1.uu();
    const int endx = nghosts + grid.dim(0);
    const int endy = nghosts + grid.dim(1);
    const int endz = nghosts + grid.dim(2);

    const double coeffx = 2. * M_PI / grid.ll(0);
    const double coeffy = 2. * M_PI / grid.ll(1);
    const double coeffz = 2. * M_PI / grid.ll(2);

    for (int ix = nghosts; ix < endx; ix++)
    {
        int iix  = ix * grid.inc(0);
        double x = grid.start(0) + ix * h[0];

        for (int iy = nghosts; iy < endy; iy++)
        {
            int iiy  = iy * grid.inc(1) + iix;
            double y = grid.start(1) + iy * h[1];

            for (int iz = nghosts; iz < endz; iz++)
            {
                double z = grid.start(2) + iz * h[2];

                u1[iiy + iz]
                    = sin(x * coeffx) + sin(y * coeffy) + sin(z * coeffz);
            }
        }
    }
    gf1.set_updated_boundaries(false);

    // fill ghost values
    gf1.trade_boundaries();

    // apply FD (-Laplacian) operator to gf1, result in gf2
    lap.apply(gf1, gf2);

    // check values in gf2
    double* u2 = gf2.uu();
    for (int ix = nghosts; ix < endx; ix++)
    {
        int iix  = ix * grid.inc(0);
        double x = grid.start(0) + ix * h[0];

        for (int iy = nghosts; iy < endy; iy++)
        {
            int iiy  = iy * grid.inc(1) + iix;
            double y = grid.start(1) + iy * h[1];

            for (int iz = nghosts; iz < endz; iz++)
            {
                double z            = grid.start(2) + iz * h[2];
                double expected_val = coeffx * coeffx * sin(x * coeffx)
                                      + coeffy * coeffy * sin(y * coeffy)
                                      + coeffz * coeffz * sin(z * coeffz);

                CHECK(u2[iiy + iz] == Approx(expected_val).margin(2.e-3));
            }
        }
    }
}
