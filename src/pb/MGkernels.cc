// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "MGkernels.h"

#include "Timer.h"
#include "mputils.h"

const double inv64 = 1. / 64.;

static Timer restrict3D_tm("restrict3D");
static Timer extend3D_tm("extend3D");

namespace pb
{

void printMGkernelTimers(std::ostream& os)
{
    restrict3D_tm.print(os);
    extend3D_tm.print(os);
}

template <typename ScalarType>
void MGkernelExtend3D(ScalarType* coarse_data, const Grid& coarse_grid,
    ScalarType* fine_data, const Grid& fine_grid, const int nfunc)
{
    extend3D_tm.start();

    assert(fine_grid.dim(0) >= 2 * coarse_grid.dim(0));
    assert(fine_grid.dim(1) >= 2 * coarse_grid.dim(1));
    assert(fine_grid.dim(2) >= 2 * coarse_grid.dim(2));
    assert(fine_grid.ghost_pt() > 0);

    const int incx_fine   = fine_grid.inc(0);
    const int incy_fine   = fine_grid.inc(1);
    const int incy_coarse = coarse_grid.inc(1);
    const int incx_coarse = coarse_grid.inc(0);

    const short nghosts_fine   = fine_grid.ghost_pt();
    const short nghosts_coarse = coarse_grid.ghost_pt();

    const int cdimx = coarse_grid.dim(0);
    const int cdimy = coarse_grid.dim(1);
    const int cdimz = coarse_grid.dim(2);

    const int dimx = fine_grid.dim(0);
    const int dimy = fine_grid.dim(1);
    const int dimz = fine_grid.dim(2);

    int ione = 1;
    int itwo = 2;

    // loop over coarse grid
    // to inject coarse values on fine grid, including first ghost
    for (int ic = 0; ic < nfunc; ic++)
    {
        int offsetf = ic * fine_grid.sizeg();
        int offsetc = ic * coarse_grid.sizeg();

        for (int ix = 0; ix <= cdimx; ix++)
        {
            const int ifx
                = offsetf + incx_fine * (2 * ix + nghosts_fine) + nghosts_fine;
            const int icx = offsetc + incx_coarse * (ix + nghosts_coarse)
                            + nghosts_coarse;

            for (int iy = 0; iy <= cdimy; iy++)
            {
                const int ify = ifx + incy_fine * (2 * iy + nghosts_fine);
                const int icy = icx + incy_coarse * (iy + nghosts_coarse);

                int size = cdimz + 1;
                Tcopy(&size, &coarse_data[icy], &ione, &fine_data[ify], &itwo);
            }
        }
    }

    // now work with fine level only
    for (int ic = 0; ic < nfunc; ic++)
    {
        int offsetf = ic * fine_grid.sizeg();

        for (int ix = nghosts_fine; ix < dimx + nghosts_fine; ix = ix + 2)
        {
            int ifx = offsetf + incx_fine * ix;

            for (int iy = nghosts_fine; iy < dimy + nghosts_fine; iy = iy + 2)
            {
                int ify = ifx + incy_fine * iy;

                for (int iz = nghosts_fine + 1; iz < dimz + nghosts_fine;
                     iz     = iz + 2)
                {
                    int izf = ify + iz;
                    fine_data[izf]
                        = 0.5 * (fine_data[izf - 1] + fine_data[izf + 1]);
                    assert(izf < static_cast<int>(nfunc * fine_grid.sizeg()));
                }
            }

            for (int iy = nghosts_fine + 1; iy < dimy + nghosts_fine;
                 iy     = iy + 2)
            {
                int ify = ifx + incy_fine * iy;

                for (int iz = nghosts_fine; iz < dimz + nghosts_fine;
                     iz     = iz + 2)
                {
                    int izf        = ify + iz;
                    fine_data[izf] = 0.5
                                     * (fine_data[izf + incy_fine]
                                         + fine_data[izf - incy_fine]);
                    assert(izf < static_cast<int>(nfunc * fine_grid.sizeg()));
                }

                for (int iz = nghosts_fine + 1; iz < dimz + nghosts_fine;
                     iz     = iz + 2)
                {
                    int izf        = ify + iz;
                    fine_data[izf] = 0.25
                                     * (fine_data[izf + 1 + incy_fine]
                                         + fine_data[izf + 1 - incy_fine]
                                         + fine_data[izf - 1 + incy_fine]
                                         + fine_data[izf - 1 - incy_fine]);
                    assert(izf < static_cast<int>(nfunc * fine_grid.sizeg()));
                }
            }
        }

        for (int ix = nghosts_fine + 1; ix < dimx + nghosts_fine; ix = ix + 2)
        {
            int ifx = offsetf + incx_fine * ix;

            for (int iy = nghosts_fine; iy < dimy + nghosts_fine; iy = iy + 2)
            {
                int ify = ifx + incy_fine * iy;

                for (int iz = nghosts_fine + 1; iz < dimz + nghosts_fine;
                     iz     = iz + 2)
                {
                    int izf        = ify + iz;
                    fine_data[izf] = 0.25
                                     * (fine_data[izf + incx_fine + 1]
                                         + fine_data[izf + incx_fine - 1]
                                         + fine_data[izf - incx_fine + 1]
                                         + fine_data[izf - incx_fine - 1]);
                }

                for (int iz = nghosts_fine; iz < dimz + nghosts_fine;
                     iz     = iz + 2)
                {
                    int izf        = ify + iz;
                    fine_data[izf] = 0.5
                                     * (fine_data[izf + incx_fine]
                                         + fine_data[izf - incx_fine]);
                }
            }

            for (int iy = nghosts_fine + 1; iy < dimy + nghosts_fine;
                 iy     = iy + 2)
            {
                int ify = ifx + incy_fine * iy;

                for (int iz = nghosts_fine + 1; iz < dimz + nghosts_fine;
                     iz     = iz + 2)
                {
                    int izf = ify + iz;
                    fine_data[izf]
                        = 0.125
                          * (fine_data[izf + incx_fine + incy_fine + 1]
                              + fine_data[izf + incx_fine + incy_fine - 1]
                              + fine_data[izf + incx_fine - incy_fine + 1]
                              + fine_data[izf + incx_fine - incy_fine - 1]
                              + fine_data[izf - incx_fine + incy_fine + 1]
                              + fine_data[izf - incx_fine + incy_fine - 1]
                              + fine_data[izf - incx_fine - incy_fine + 1]
                              + fine_data[izf - incx_fine - incy_fine - 1]);
                }

                for (int iz = nghosts_fine; iz < dimz + nghosts_fine;
                     iz     = iz + 2)
                {
                    int izf = ify + iz;
                    fine_data[izf]
                        = 0.25
                          * (fine_data[izf + incx_fine + incy_fine]
                              + fine_data[izf + incx_fine - incy_fine]
                              + fine_data[izf - incx_fine + incy_fine]
                              + fine_data[izf - incx_fine - incy_fine]);
                }
            }
        }
    }

    extend3D_tm.stop();
}

template <typename ScalarType>
void MGkernelRestrict3D(ScalarType* fine_data, const Grid& fine_grid,
    ScalarType* coarse_data, const Grid& coarse_grid, const int nfunc)
{
    restrict3D_tm.start();

    const int incx_fine   = fine_grid.inc(0);
    const int incy_fine   = fine_grid.inc(1);
    const int incy_coarse = coarse_grid.inc(1);
    const int incx_coarse = coarse_grid.inc(0);

    const short nghosts_fine   = fine_grid.ghost_pt();
    const short nghosts_coarse = coarse_grid.ghost_pt();

    const int cdimx = coarse_grid.dim(0);
    const int cdimy = coarse_grid.dim(1);
    const int cdimz = coarse_grid.dim(2);

    // loop over coarse grid
    for (int ifunc = 0; ifunc < nfunc; ifunc++)
    {
        int ixf = ifunc * fine_grid.sizeg() + nghosts_fine * incx_fine;
        int ixc = ifunc * coarse_grid.sizeg() + nghosts_coarse * incx_coarse;

        for (int ix = nghosts_fine; ix < cdimx + nghosts_fine; ix++)
        {
            int iyf = ixf + incy_fine * nghosts_fine;
            int iyc = ixc + incy_coarse * nghosts_coarse;

            for (int iy = nghosts_fine; iy < cdimy + nghosts_fine; iy++)
            {
                int izf = iyf + nghosts_fine;

                const ScalarType* const u0    = fine_data + izf;
                const ScalarType* const umx   = u0 - incx_fine;
                const ScalarType* const upx   = u0 + incx_fine;
                const ScalarType* const umy   = u0 - incy_fine;
                const ScalarType* const upy   = u0 + incy_fine;
                const ScalarType* const umxpy = u0 - incx_fine + incy_fine;
                const ScalarType* const upxpy = u0 + incx_fine + incy_fine;
                const ScalarType* const umxmy = u0 - incx_fine - incy_fine;
                const ScalarType* const upxmy = u0 + incx_fine - incy_fine;
                for (int iz = 0; iz < cdimz; iz++)
                {
                    const int twoiz = 2 * iz;

                    double face = (double)upx[twoiz] + (double)umx[twoiz]
                                  + (double)upy[twoiz] + (double)umy[twoiz]
                                  + (double)u0[twoiz - 1]
                                  + (double)u0[twoiz + 1];

                    double corner
                        = (double)upxpy[twoiz - 1] + (double)upxpy[twoiz + 1]
                          + (double)upxmy[twoiz - 1] + (double)upxmy[twoiz + 1]
                          + (double)umxpy[twoiz - 1] + (double)umxpy[twoiz + 1]
                          + (double)umxmy[twoiz - 1] + (double)umxmy[twoiz + 1];

                    double edge
                        = (double)upy[twoiz - 1] + (double)upy[twoiz + 1]
                          + (double)umy[twoiz - 1] + (double)umy[twoiz + 1]
                          + (double)upx[twoiz - 1] + (double)upx[twoiz + 1]
                          + (double)umx[twoiz - 1] + (double)umx[twoiz + 1]
                          + (double)umxmy[twoiz] + (double)upxmy[twoiz]
                          + (double)umxpy[twoiz] + (double)upxpy[twoiz];

                    coarse_data[iyc + iz + nghosts_coarse]
                        = (ScalarType)(inv64
                                       * (8. * u0[twoiz] + 4. * face + 2. * edge
                                           + corner));
                }

                iyf += 2 * incy_fine;
                iyc += incy_coarse;
            }

            ixf += 2 * incx_fine;
            ixc += incx_coarse;
        }
    }

    restrict3D_tm.stop();
}

template void MGkernelExtend3D<double>(double* coarse_data,
    const Grid& coarse_grid, double* fine_data, const Grid& fine_grid,
    const int nfunc);
template void MGkernelExtend3D<float>(float* coarse_data,
    const Grid& coarse_grid, float* fine_data, const Grid& fine_grid,
    const int nfunc);
template void MGkernelRestrict3D<double>(double* fine_data,
    const Grid& fine_grid, double* coarse_data, const Grid& coarse_grid,
    const int nfunc);
template void MGkernelRestrict3D<float>(float* fine_data, const Grid& fine_grid,
    float* coarse_data, const Grid& coarse_grid, const int nfunc);
}
