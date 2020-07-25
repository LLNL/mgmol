#include "GridFunc.h"
#include "Laph4.h"
#include "PEenv.h"

#include "catch.hpp"

TEST_CASE("Laplacian 4th order for a array of function", "[laph4_batch]")
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

    const size_t numgridfunc = 10;

    double* arrayofgf1 = new double[numgridfunc * grid.sizeg()];

    const int endx = nghosts + grid.dim(0);
    const int endy = nghosts + grid.dim(1);
    const int endz = nghosts + grid.dim(2);

    const double coeffx = 2. * M_PI / grid.ll(0);
    const double coeffy = 2. * M_PI / grid.ll(1);
    const double coeffz = 2. * M_PI / grid.ll(2);

    // periodic GridFunc
    pb::GridFunc<double> gf1(grid, 1, 1, 1);
    pb::GridFunc<double> gf2(grid, 1, 1, 1);

    double* u1 = gf1.uu();

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

    for (size_t inum = 0; inum < numgridfunc; inum++)
    {
        std::copy(gf1.uu(), gf1.uu() + grid.sizeg(),
            arrayofgf1 + inum * grid.sizeg());
    }

    double* arrayofgf2 = new double[numgridfunc * grid.sizeg()];

    // apply FD (-Laplacian) operator to arrayofgf1, result in arrayofgf2
    lap.apply(grid, arrayofgf1, arrayofgf2, numgridfunc);

    // check values in gf2
    double* u2 = gf2.uu();

    for (size_t inum = 0; inum < numgridfunc; inum++)
    {
        std::copy(arrayofgf2 + inum * grid.sizeg(),
            arrayofgf2 + (inum + 1) * grid.sizeg(), u2);

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

    delete arrayofgf1;
    delete arrayofgf2;
}
