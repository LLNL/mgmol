#include "FDkernels.h"
#include "Timer.h"

static Timer del2_2nd_tm("del2_2nd");

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

template void FDkernelDel2_2nd<double>(const Grid& grid, double* v, double* b,
    const size_t nfunc, MemorySpace::Host);
template void FDkernelDel2_2nd<float>(const Grid& grid, float* v, float* b,
    const size_t nfunc, MemorySpace::Host);

} // namespace pb
