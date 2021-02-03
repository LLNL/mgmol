#include "FDkernels.h"
#include "Timer.h"

static Timer del2_2nd_tm("del2_2nd");
static Timer del2_4th_tm("del2_4th");

const double inv12 = 1. / 12.;

namespace pb
{

void printFDkernelTimers(std::ostream& os) { del2_2nd_tm.print(os); }

// evaluate 2nd order FD Laplacian on uniform mesh defined by Grid
// object
template <typename ScalarType>
void FDkernelDel2_2nd(const Grid& grid, ScalarType* v, ScalarType* b,
    const size_t nfunc, MemorySpace::Host)
{
    assert(grid.ghost_pt() > 0);
    MemorySpace::assert_is_host_ptr(b);
    MemorySpace::assert_is_host_ptr(v);

    del2_2nd_tm.start();

    double inv_h2[3] = {
        1. / (grid.hgrid(0) * grid.hgrid(0)),
        1. / (grid.hgrid(1) * grid.hgrid(1)),
        1. / (grid.hgrid(2) * grid.hgrid(2)),
    };

    double cc        = inv_h2[0];
    const double c1x = -cc;

    cc               = inv_h2[1];
    const double c1y = -cc;

    cc               = inv_h2[2];
    const double c1z = -cc;
    const double c0  = -2. * (c1x + c1y + c1z);

    const int dim0 = grid.dim(0);
    const int dim1 = grid.dim(1);
    const int dim2 = grid.dim(2);

    const size_t ngpts = grid.sizeg();

    int incx     = grid.inc(0);
    int incy     = grid.inc(1);
    const int ng = grid.ghost_pt();

#pragma omp parallel for collapse(4)
    for (size_t ifunc = 0; ifunc < nfunc; ifunc++)
    {
        for (int ix = 0; ix < dim0; ix++)
        {
            for (int iy = 0; iy < dim1; iy++)
            {
                for (int iz = 0; iz < dim2; iz++)
                {
                    int iiz = ifunc * ngpts
                              + ((ix + ng) * incx + (iy + ng) * incy) + iz + ng;

                    b[iiz] = c0 * v[iiz] + c1x * (v[iiz - incx] + v[iiz + incx])
                             + c1y * (v[iiz - incy] + v[iiz + incy])
                             + c1z * (v[iiz - 1] + v[iiz + 1]);
                }
            }
        }
    }

    del2_2nd_tm.stop();
}

template <typename ScalarType>
void FDkernelDel2_4th(const Grid& grid, ScalarType* A, ScalarType* B,
    const size_t nfunc, MemorySpace::Host)
{
    MemorySpace::assert_is_host_ptr(A);
    MemorySpace::assert_is_host_ptr(B);

    del2_4th_tm.start();

    double inv_h2[3] = {
        1. / (grid.hgrid(0) * grid.hgrid(0)),
        1. / (grid.hgrid(1) * grid.hgrid(1)),
        1. / (grid.hgrid(2) * grid.hgrid(2)),
    };

    const double cc0 = inv12 * inv_h2[0];
    const double c1x = -16. * cc0;
    const double c2x = 1. * cc0;

    const double cc1 = inv12 * inv_h2[1];
    const double c1y = -16. * cc1;
    const double c2y = 1. * cc1;

    const double cc2 = inv12 * inv_h2[2];
    const double c1z = -16. * cc2;
    const double c2z = 1. * cc2;

    const double c0 = -2. * (c1x + c2x + c1y + c2y + c1z + c2z);

    const int gpt = grid.ghost_pt();

    const int incx = grid.inc(0);
    const int incy = grid.inc(1);

    const int incx2 = 2 * grid.inc(0);
    const int incy2 = 2 * grid.inc(1);

    int iix = gpt * incx;

    const int dim0 = grid.dim(0);
    const int dim1 = grid.dim(1);
    const int dim2 = grid.dim(2);

    const size_t ng = grid.sizeg();

#pragma omp parallel for collapse(4)
    for (size_t ifunc = 0; ifunc < nfunc; ifunc++)
    {
        for (int ix = 0; ix < dim0; ix++)
        {
            for (int iy = 0; iy < dim1; iy++)
            {
                for (int iz = 0; iz < dim2; iz++)
                {
                    int iiz = ifunc * ng
                              + (iix + ix * incx + gpt * incy + iy * incy)
                              + gpt;

                    const ScalarType* __restrict__ v0   = A + iiz;
                    const ScalarType* __restrict__ vmx  = A + (iiz - incx);
                    const ScalarType* __restrict__ vpx  = A + (iiz + incx);
                    const ScalarType* __restrict__ vmx2 = A + (iiz - incx2);
                    const ScalarType* __restrict__ vpx2 = A + (iiz + incx2);
                    const ScalarType* __restrict__ vmy  = A + (iiz - incy);
                    const ScalarType* __restrict__ vpy  = A + (iiz + incy);
                    const ScalarType* __restrict__ vmy2 = A + (iiz - incy2);
                    const ScalarType* __restrict__ vpy2 = A + (iiz + incy2);

                    ScalarType* __restrict__ u = B + iiz;

                    u[iz] = c0 * (double)v0[iz]

                            + c1x * ((double)vmx[iz] + (double)vpx[iz])
                            + c1y * ((double)vmy[iz] + (double)vpy[iz])
                            + c1z * ((double)v0[iz - 1] + (double)v0[iz + 1])

                            + c2x * ((double)vmx2[iz] + (double)vpx2[iz])
                            + c2y * ((double)vmy2[iz] + (double)vpy2[iz])
                            + c2z * ((double)v0[iz - 2] + (double)v0[iz + 2]);
                }
            }
        }
    }

    del2_4th_tm.stop();
}

#ifdef HAVE_MAGMA
template <typename ScalarType>
void FDkernelDel2_4th(const Grid& grid, ScalarType* A, ScalarType* B,
    const size_t nfunc, MemorySpace::Device)
{
    MemorySpace::assert_is_dev_ptr(A);
    MemorySpace::assert_is_dev_ptr(B);

    assert(grid.ghost_pt() > 1);

    del2_4th_tm.start();

    double inv_h2[3] = {
        1. / (grid.hgrid(0) * grid.hgrid(0)),
        1. / (grid.hgrid(1) * grid.hgrid(1)),
        1. / (grid.hgrid(2) * grid.hgrid(2)),
    };

    const double cc0 = inv12 * inv_h2[0];
    const double c1x = -16. * cc0;
    const double c2x = 1. * cc0;

    const double cc1 = inv12 * inv_h2[1];
    const double c1y = -16. * cc1;
    const double c2y = 1. * cc1;

    const double cc2 = inv12 * inv_h2[2];
    const double c1z = -16. * cc2;
    const double c2z = 1. * cc2;

    const double c0 = -2. * (c1x + c2x + c1y + c2y + c1z + c2z);

    const int gpt = grid.ghost_pt();

    const int incx = grid.inc(0);
    const int incy = grid.inc(1);

    const int incx2 = 2 * grid.inc(0);
    const int incy2 = 2 * grid.inc(1);

    int iix = gpt * incx_;

    const int dim0 = Agrid.dim(0);
    const int dim1 = Agrid.dim(1);
    const int dim2 = Agrid.dim(2);

    const size_t ng = grid.sizeg();

    MGMOL_PARALLEL_FOR_COLLAPSE(2, A, B) // 0.037 s 3000x32^3 0.0156 1000x32^3
    for (size_t ifunc = 0; ifunc < nfunc; ifunc++)
    {
        for (int ix = 0; ix < dim0; ix++)
        {
            for (int iy = 0; iy < dim1; iy++)
            {
                for (int iz = 0; iz < dim2; iz++)
                {
                    int iiz = ifunc * ng
                              + (iix + ix * incx + gpt * incy + iy * incy)
                              + gpt;

                    const ScalarType* __restrict__ v0   = A + iiz;
                    const ScalarType* __restrict__ vmx  = A + (iiz - incx);
                    const ScalarType* __restrict__ vpx  = A + (iiz + incx);
                    const ScalarType* __restrict__ vmx2 = A + (iiz - incx2);
                    const ScalarType* __restrict__ vpx2 = A + (iiz + incx2);
                    const ScalarType* __restrict__ vmy  = A + (iiz - incy);
                    const ScalarType* __restrict__ vpy  = A + (iiz + incy);
                    const ScalarType* __restrict__ vmy2 = A + (iiz - incy2);
                    const ScalarType* __restrict__ vpy2 = A + (iiz + incy2);

                    ScalarType* __restrict__ u = B + iiz;

                    u[iz] = c0 * (double)v0[iz]

                            + c1x * ((double)vmx[iz] + (double)vpx[iz])
                            + c1y * ((double)vmy[iz] + (double)vpy[iz])
                            + c1z * ((double)v0[iz - 1] + (double)v0[iz + 1])

                            + c2x * ((double)vmx2[iz] + (double)vpx2[iz])
                            + c2y * ((double)vmy2[iz] + (double)vpy2[iz])
                            + c2z * ((double)v0[iz - 2] + (double)v0[iz + 2]);
                }
            }
        }
    }

    del2_4th_tm.stop();
}
#endif

template void FDkernelDel2_2nd<double>(const Grid& grid, double* v, double* b,
    const size_t nfunc, MemorySpace::Host);
template void FDkernelDel2_2nd<float>(const Grid& grid, float* v, float* b,
    const size_t nfunc, MemorySpace::Host);

template void FDkernelDel2_4th<double>(const Grid& grid, double* v, double* b,
    const size_t nfunc, MemorySpace::Host);
template void FDkernelDel2_4th<float>(const Grid& grid, float* v, float* b,
    const size_t nfunc, MemorySpace::Host);

#ifdef HAVE_MAGMA
template void FDkernelDel2_4th<double>(const Grid& grid, double* v, double* b,
    const size_t nfunc, MemorySpace::Device);
template void FDkernelDel2_4th<float>(const Grid& grid, float* v, float* b,
    const size_t nfunc, MemorySpace::Device);
#endif
} // namespace pb
