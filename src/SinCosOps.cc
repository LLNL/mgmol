// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SinCosOps.h"

#include "ExtendedGridOrbitals.h"
#include "FunctionsPacking.h"
#include "LocGridOrbitals.h"
#include "MGmol_MPI.h"
#include "memory_space.h"

using namespace std;

#ifdef HAVE_MAGMA
using memory_space_type = MemorySpace::Device;
#else
using memory_space_type = MemorySpace::Host;
#endif

template <class T>
void SinCosOps<T>::compute(const T& orbitals, vector<vector<double>>& a)
{
    assert(a.size() == 6);

    compute_tm_.start();

    const pb::Grid& grid(orbitals.grid_);
    const int numst = orbitals.numst();

    const int dim0 = grid.dim(0);
    const int dim1 = grid.dim(1);
    const int dim2 = grid.dim(2);

    int n2 = numst * numst;

    int loc_length = dim0 / orbitals.subdivx();
    assert(loc_length > 0);
    assert(loc_length <= dim0);

    int incx = dim1 * dim2;
    int incy = dim2;

    vector<double> sinx;
    vector<double> siny;
    vector<double> sinz;
    vector<double> cosx;
    vector<double> cosy;
    vector<double> cosz;
    grid.getSinCosFunctions(sinx, siny, sinz, cosx, cosy, cosz);

    const int size = orbitals.chromatic_number();

    const int ld                = orbitals.getLda();
    unsigned int const size_psi = size * ld;
    ORBDTYPE* psi_view
        = MemorySpace::Memory<ORBDTYPE, memory_space_type>::allocate_host_view(
            size_psi);
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
        orbitals.psi(0), size_psi, psi_view);

    for (short iloc = 0; iloc < orbitals.subdivx(); iloc++)
    {

        for (int icolor = 0; icolor < size; icolor++)
        {
            const int i = orbitals.overlapping_gids_[iloc][icolor];
            if (i != -1)
            {
                const ORBDTYPE* const ppsii = psi_view + ld * icolor;
                for (int jstate = 0; jstate <= icolor; jstate++)
                {
                    const int j = orbitals.overlapping_gids_[iloc][jstate];
                    if (j != -1)
                    {

                        const ORBDTYPE* const ppsij = psi_view + ld * jstate;

                        double atmp[6]  = { 0., 0., 0., 0., 0., 0. };
                        const int ixend = loc_length * (iloc + 1);

                        for (int ix = loc_length * iloc; ix < ixend; ix++)
                        {
                            const double cosix = cosx[ix];
                            const double sinix = sinx[ix];
                            const int offsetx  = ix * incx;

                            for (int iy = 0; iy < dim1; iy++)
                            {
                                const double cosiy = cosy[iy];
                                const double siniy = siny[iy];
                                const int offset   = offsetx + iy * incy;

                                for (int iz = 0; iz < dim2; iz++)
                                {
                                    const int index    = offset + iz;
                                    const double alpha = (double)ppsij[index]
                                                         * (double)ppsii[index];
                                    atmp[0] += alpha * cosix;
                                    atmp[1] += alpha * sinix;
                                    atmp[2] += alpha * cosiy;
                                    atmp[3] += alpha * siniy;
                                    atmp[4] += alpha * cosz[iz];
                                    atmp[5] += alpha * sinz[iz];
                                }
                            }
                        }
                        const int ji = j * numst + i;
                        const int ij = i * numst + j;
                        a[0][ji]     = a[0][ij] += atmp[0];
                        a[1][ji]     = a[1][ij] += atmp[1];
                        a[2][ji]     = a[2][ij] += atmp[2];
                        a[3][ji]     = a[3][ij] += atmp[3];
                        a[4][ji]     = a[4][ij] += atmp[4];
                        a[5][ji]     = a[5][ij] += atmp[5];
                    }
                }
            }
        }
    }

    MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(psi_view);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    for (short i = 0; i < 6; i++)
    {
        mmpi.split_allreduce_sums_double(&a[i][0], n2);
        my_dscal(n2, grid.vel(), &a[i][0]);
    }

    compute_tm_.stop();
}

template <class T>
void SinCosOps<T>::computeSquare(const T& orbitals, vector<vector<double>>& a)
{
    assert(a.size() == 6);
    for (short i = 0; i < 6; i++)
        assert(a[i].size() > 0);

    compute_tm_.start();

    const pb::Grid& grid(orbitals.grid_);
    const int numst = orbitals.numst_;
    const int n2    = numst * numst;

    const int dim0 = grid.dim(0);
    const int dim1 = grid.dim(1);
    const int dim2 = grid.dim(2);

    int loc_length = dim0 / orbitals.subdivx();
    assert(loc_length > 0);
    assert(loc_length <= dim0);

    const int xoff = grid.istart(0);
    const int yoff = grid.istart(1);
    const int zoff = grid.istart(2);

    const double hhx = 2. * M_PI / (double)grid.gdim(0);
    const double hhy = 2. * M_PI / (double)grid.gdim(1);
    const double hhz = 2. * M_PI / (double)grid.gdim(2);

    int incx = dim1 * dim2;
    int incy = dim2;

    const double inv_2pi = 0.5 * M_1_PI;
    const double alphax  = grid.ll(0) * grid.ll(0) * inv_2pi * inv_2pi;
    const double alphay  = grid.ll(1) * grid.ll(1) * inv_2pi * inv_2pi;
    const double alphaz  = grid.ll(2) * grid.ll(2) * inv_2pi * inv_2pi;

    vector<double> sinx2, siny2, sinz2, cosx2, cosy2, cosz2;
    sinx2.resize(dim0);
    cosx2.resize(dim0);
    siny2.resize(dim1);
    cosy2.resize(dim1);
    sinz2.resize(dim2);
    cosz2.resize(dim2);
    for (int i = 0; i < dim0; i++)
    {
        const double tmp = sin((double)(xoff + i) * hhx);
        sinx2[i]         = tmp * tmp * alphax;
        cosx2[i]         = (1. - tmp * tmp) * alphax;
    }
    for (int i = 0; i < dim1; i++)
    {
        const double tmp = sin((double)(yoff + i) * hhy);
        siny2[i]         = tmp * tmp * alphay;
        cosy2[i]         = (1. - tmp * tmp) * alphay;
    }
    for (int i = 0; i < dim2; i++)
    {
        const double tmp = sin((double)(zoff + i) * hhz);
        sinz2[i]         = tmp * tmp * alphaz;
        cosz2[i]         = (1. - tmp * tmp) * alphaz;
    }
    const int size = orbitals.chromatic_number();

    for (short iloc = 0; iloc < orbitals.subdivx(); iloc++)
    {

        for (int icolor = 0; icolor < size; icolor++)
        {
            int i = orbitals.overlapping_gids_[iloc][icolor];
            if (i != -1)
                for (int jstate = 0; jstate <= icolor; jstate++)
                {
                    int j = orbitals.overlapping_gids_[iloc][jstate];
                    if (j != -1)
                    {

                        double atmp[6] = { 0., 0., 0., 0., 0., 0. };
                        for (int ix = loc_length * iloc;
                             ix < loc_length * (iloc + 1); ix++)
                            for (int iy = 0; iy < dim1; iy++)
                                for (int iz = 0; iz < dim2; iz++)
                                {

                                    int index = ix * incx + iy * incy + iz;
                                    double alpha
                                        = orbitals.psi(jstate)[index]
                                          * orbitals.psi(icolor)[index];
                                    atmp[0] += alpha * cosx2[ix];
                                    atmp[1] += alpha * sinx2[ix];
                                    atmp[2] += alpha * cosy2[iy];
                                    atmp[3] += alpha * siny2[iy];
                                    atmp[4] += alpha * cosz2[iz];
                                    atmp[5] += alpha * sinz2[iz];
                                }
                        int ji   = j * numst + i;
                        int ij   = i * numst + j;
                        a[0][ji] = a[0][ij] += atmp[0];
                        a[1][ji] = a[1][ij] += atmp[1];
                        a[2][ji] = a[2][ij] += atmp[2];
                        a[3][ji] = a[3][ij] += atmp[3];
                        a[4][ji] = a[4][ij] += atmp[4];
                        a[5][ji] = a[5][ij] += atmp[5];
                    }
                }
        }
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    for (short i = 0; i < 6; i++)
    {
        mmpi.split_allreduce_sums_double(&a[i][0], n2);
        my_dscal(n2, grid.vel(), &a[i][0]);
    }

    compute_tm_.stop();
}

template <class T>
void SinCosOps<T>::computeSquare1D(
    const T& orbitals, vector<vector<double>>& a, const int dim_index)
{
    assert(a.size() == 2);
    for (short i = 0; i < 2; i++)
        assert(a[i].size() > 0);

    compute_tm_.start();

    const pb::Grid& grid(orbitals.grid_);
    const int numst = orbitals.numst_;

    const int dim  = grid.dim(dim_index);
    const int dim0 = grid.dim(0);
    const int dim1 = grid.dim(1);
    const int dim2 = grid.dim(2);

    int n2 = numst * numst;

    int loc_length = dim0 / orbitals.subdivx();
    assert(loc_length > 0);
    assert(loc_length <= dim0);

    const int off   = grid.istart(dim_index);
    const double hh = 2. * M_PI / (double)grid.gdim(dim_index);

    int incx = dim1 * dim2;
    int incy = dim2;

    const double inv_2pi = 0.5 * M_1_PI;
    const double alphax
        = grid.ll(dim_index) * grid.ll(dim_index) * inv_2pi * inv_2pi;

    vector<double> sinx2(dim);
    vector<double> cosx2(dim);
    for (int i = 0; i < dim; i++)
    {
        const double tmp = sin((double)(off + i) * hh);
        sinx2[i]         = tmp * tmp * alphax;
        cosx2[i]         = (1. - tmp * tmp) * alphax;
    }
    const int size = orbitals.chromatic_number();

    for (short iloc = 0; iloc < orbitals.subdivx(); iloc++)
    {

        for (int icolor = 0; icolor < size; icolor++)
        {
            int i = orbitals.overlapping_gids_[iloc][icolor];
            if (i != -1)
                for (int jstate = 0; jstate <= icolor; jstate++)
                {
                    int j = orbitals.overlapping_gids_[iloc][jstate];
                    if (j != -1)
                    {

                        double atmp[2] = { 0., 0. };
                        for (int ix = loc_length * iloc;
                             ix < loc_length * (iloc + 1); ix++)
                            for (int iy = 0; iy < dim1; iy++)
                                for (int iz = 0; iz < dim2; iz++)
                                {

                                    int dindex[3] = { ix, iy, iz };
                                    int index     = ix * incx + iy * incy + iz;
                                    double alpha
                                        = orbitals.psi(jstate)[index]
                                          * orbitals.psi(icolor)[index];
                                    atmp[0] += alpha * cosx2[dindex[dim_index]];
                                    atmp[1] += alpha * sinx2[dindex[dim_index]];
                                }
                        int ji   = j * numst + i;
                        int ij   = i * numst + j;
                        a[0][ji] = a[0][ij] += atmp[0];
                        a[1][ji] = a[1][ij] += atmp[1];
                    }
                }
        }
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    for (int i = 0; i < 2; i++)
    {
        mmpi.split_allreduce_sums_double(&a[i][0], n2);
        my_dscal(n2, grid.vel(), &a[i][0]);
    }

    compute_tm_.stop();
}

template <class T>
void SinCosOps<T>::compute1D(
    const T& orbitals, vector<vector<double>>& a, const int dim_index)
{
    assert(a.size() == 2);
    for (short i = 0; i < 2; i++)
        assert(a[i].size() > 0);

    compute_tm_.start();

    const pb::Grid& grid(orbitals.grid_);
    const int numst = orbitals.numst_;

    const int dim = grid.dim(dim_index);

    const int dim0 = grid.dim(0);
    const int dim1 = grid.dim(1);
    const int dim2 = grid.dim(2);

    int n2 = numst * numst;

    int loc_length = dim0 / orbitals.subdivx();
    assert(loc_length > 0);
    assert(loc_length <= dim0);

    int incx = dim1 * dim2;
    int incy = dim2;

    const int off   = grid.istart(dim_index);
    const double hh = 2. * M_PI / (double)grid.gdim(dim_index);

    const double inv_2pi = 0.5 * M_1_PI;
    const double alphax  = grid.ll(dim_index) * inv_2pi;

    vector<double> sinx(dim);
    for (int i = 0; i < dim; i++)
        sinx[i] = sin(double(off + i) * hh) * alphax;

    vector<double> cosx(dim);
    for (int i = 0; i < dim; i++)
        cosx[i] = cos(double(off + i) * hh) * alphax;

    const int size = orbitals.chromatic_number();

    for (short iloc = 0; iloc < orbitals.subdivx(); iloc++)
    {

        for (int icolor = 0; icolor < size; icolor++)
        {
            const int i = orbitals.overlapping_gids_[iloc][icolor];
            if (i != -1)
            {
                const ORBDTYPE* const ppsii = orbitals.psi(icolor);
                for (int jstate = 0; jstate <= icolor; jstate++)
                {
                    const int j = orbitals.overlapping_gids_[iloc][jstate];
                    if (j != -1)
                    {

                        const ORBDTYPE* const ppsij = orbitals.psi(jstate);

                        double atmp[2]  = { 0., 0. };
                        const int ixend = loc_length * (iloc + 1);

                        for (int ix = loc_length * iloc; ix < ixend; ix++)
                            for (int iy = 0; iy < dim1; iy++)
                                for (int iz = 0; iz < dim2; iz++)
                                {

                                    int dindex[3] = { ix, iy, iz };
                                    const int index
                                        = ix * incx + iy * incy + iz;
                                    const double alpha = (double)ppsij[index]
                                                         * (double)ppsii[index];
                                    atmp[0] += alpha * cosx[dindex[dim_index]];
                                    atmp[1] += alpha * sinx[dindex[dim_index]];
                                }
                        const int ji = j * numst + i;
                        const int ij = i * numst + j;
                        a[0][ji]     = a[0][ij] += atmp[0];
                        a[1][ji]     = a[1][ij] += atmp[1];
                    }
                }
            }
        }
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    for (short i = 0; i < 2; i++)
    {
        mmpi.split_allreduce_sums_double(&a[i][0], n2);
        my_dscal(n2, grid.vel(), &a[i][0]);
    }

    compute_tm_.stop();
}

template <class T>
void SinCosOps<T>::computeDiag2states(
    const T& orbitals, vector<vector<double>>& a, const int st1, const int st2)
{
    assert(st1 >= 0);
    assert(st2 >= 0);
    assert(st1 != st2);

    compute_tm_.start();

    const pb::Grid& grid(orbitals.grid_);

    const int st[2] = { st1, st2 };

    const int dim0 = grid.dim(0);
    const int dim1 = grid.dim(1);
    const int dim2 = grid.dim(2);

    short color_st[2] = { -1, -1 };
    for (short ic = 0; ic < 2; ++ic)
    {
        color_st[ic] = orbitals.getColor(st[ic]);
    }

    int loc_length = dim0 / orbitals.subdivx();
    assert(loc_length > 0);
    assert(loc_length <= dim0);

    int incx = dim1 * dim2;
    int incy = dim2;

    vector<double> sinx;
    vector<double> siny;
    vector<double> sinz;
    vector<double> cosx;
    vector<double> cosy;
    vector<double> cosz;
    grid.getSinCosFunctions(sinx, siny, sinz, cosx, cosy, cosz);

    double norm2[2] = { 0., 0. };
    for (short ic = 0; ic < 2; ic++)
    {
        const short mycolor = color_st[ic];

        if (mycolor >= 0)
            for (short iloc = 0; iloc < orbitals.subdivx(); iloc++)
            {

                if (orbitals.overlapping_gids_[iloc][mycolor] == st[ic])
                {
                    const ORBDTYPE* const ppsii = orbitals.psi(mycolor);
                    assert(ppsii != nullptr);
                    double atmp[6]  = { 0., 0., 0., 0., 0., 0. };
                    const int ixend = loc_length * (iloc + 1);

                    for (int ix = loc_length * iloc; ix < ixend; ix++)
                        for (int iy = 0; iy < dim1; iy++)
                            for (int iz = 0; iz < dim2; iz++)
                            {

                                const int index    = ix * incx + iy * incy + iz;
                                const double alpha = (double)ppsii[index]
                                                     * (double)ppsii[index];
                                atmp[0] += alpha * cosx[ix];
                                atmp[1] += alpha * sinx[ix];
                                atmp[2] += alpha * cosy[iy];
                                atmp[3] += alpha * siny[iy];
                                atmp[4] += alpha * cosz[iz];
                                atmp[5] += alpha * sinz[iz];
                                norm2[ic] += alpha;
                            }
                    a[0][ic] += atmp[0];
                    a[1][ic] += atmp[1];
                    a[2][ic] += atmp[2];
                    a[3][ic] += atmp[3];
                    a[4][ic] += atmp[4];
                    a[5][ic] += atmp[5];
                }
            }
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    for (short i = 0; i < 6; i++)
    {
        mmpi.split_allreduce_sums_double(&a[i][0], 2);
        my_dscal(2, grid.vel(), &a[i][0]);
    }

    compute_tm_.stop();
}

template <class T>
void SinCosOps<T>::compute2states(
    const T& orbitals, vector<vector<double>>& a, const int st1, const int st2)
{
    assert(a.size() == 6);
    assert(st1 >= 0);
    assert(st2 >= 0);

    compute_tm_.start();

    const pb::Grid& grid(orbitals.grid_);

    const int st[2] = { st1, st2 };
    int color_st[2] = { -1, -1 };
    for (short ic = 0; ic < 2; ++ic)
    {
        color_st[ic] = orbitals.getColor(st[ic]);
    }

    const int dim0 = grid.dim(0);
    const int dim1 = grid.dim(1);
    const int dim2 = grid.dim(2);

    int n2 = 4;

    int loc_length = dim0 / orbitals.subdivx();
    assert(loc_length > 0);
    assert(loc_length <= dim0);

    int incx = dim1 * dim2;
    int incy = dim2;

    vector<double> sinx;
    vector<double> siny;
    vector<double> sinz;
    vector<double> cosx;
    vector<double> cosy;
    vector<double> cosz;
    grid.getSinCosFunctions(sinx, siny, sinz, cosx, cosy, cosz);

    for (int ic = 0; ic < 2; ic++)
    {
        const int mycolor = color_st[ic];

        if (mycolor >= 0)
            for (short iloc = 0; iloc < orbitals.subdivx(); iloc++)
            {

                if (orbitals.overlapping_gids_[iloc][mycolor] == st[ic])
                {
                    const ORBDTYPE* const ppsii = orbitals.psi(mycolor);
                    assert(ppsii != nullptr);
                    for (int jc = 0; jc <= ic; jc++)
                        if (color_st[jc] >= 0)
                        {
                            if (orbitals.overlapping_gids_[iloc][color_st[jc]]
                                == st[jc])
                            {
                                const ORBDTYPE* const ppsij
                                    = orbitals.psi(color_st[jc]);
                                assert(ppsij != nullptr);

                                double atmp[6]  = { 0., 0., 0., 0., 0., 0. };
                                const int ixend = loc_length * (iloc + 1);

                                for (int ix = loc_length * iloc; ix < ixend;
                                     ix++)
                                    for (int iy = 0; iy < dim1; iy++)
                                        for (int iz = 0; iz < dim2; iz++)
                                        {

                                            const int index
                                                = ix * incx + iy * incy + iz;
                                            const double alpha
                                                = (double)ppsij[index]
                                                  * (double)ppsii[index];
                                            atmp[0] += alpha * cosx[ix];
                                            atmp[1] += alpha * sinx[ix];
                                            atmp[2] += alpha * cosy[iy];
                                            atmp[3] += alpha * siny[iy];
                                            atmp[4] += alpha * cosz[iz];
                                            atmp[5] += alpha * sinz[iz];
                                        }
                                const int ji = jc * 2 + ic;
                                const int ij = ic * 2 + jc;
                                a[0][ji]     = a[0][ij] += atmp[0];
                                a[1][ji]     = a[1][ij] += atmp[1];
                                a[2][ji]     = a[2][ij] += atmp[2];
                                a[3][ji]     = a[3][ij] += atmp[3];
                                a[4][ji]     = a[4][ij] += atmp[4];
                                a[5][ji]     = a[5][ij] += atmp[5];
                            }
                        }
                }
            }
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    for (short i = 0; i < 6; i++)
    {
        mmpi.split_allreduce_sums_double(&a[i][0], n2);
        my_dscal(n2, grid.vel(), &a[i][0]);
    }

    compute_tm_.stop();
}

template <class T>
void SinCosOps<T>::compute(
    const T& orbitals1, const T& orbitals2, vector<vector<double>>& a)
{
    assert(a.size() == 6);

    compute_tm_.start();

    const pb::Grid& grid(orbitals1.grid_);
    const int numst = orbitals1.numst_;
    int n2          = numst * numst;

    const int dim0 = grid.dim(0);
    const int dim1 = grid.dim(1);
    const int dim2 = grid.dim(2);

    int loc_length = dim0 / orbitals1.subdivx();
    assert(loc_length > 0);
    assert(loc_length <= dim0);

    int incx = dim1 * dim2;
    int incy = dim2;

    vector<double> sinx;
    vector<double> siny;
    vector<double> sinz;
    vector<double> cosx;
    vector<double> cosy;
    vector<double> cosz;
    grid.getSinCosFunctions(sinx, siny, sinz, cosx, cosy, cosz);

    for (short iloc = 0; iloc < orbitals1.subdivx(); iloc++)
    {
        for (int color = 0; color < orbitals1.chromatic_number(); color++)
        {
            int i = orbitals1.overlapping_gids_[iloc][color];
            if (i != -1)
                for (int jstate = 0; jstate < orbitals2.chromatic_number();
                     jstate++)
                {
                    int j = orbitals2.overlapping_gids_[iloc][jstate];
                    if (j != -1)
                    {

                        double atmp[6] = { 0., 0., 0., 0., 0., 0. };
                        for (int ix = loc_length * iloc;
                             ix < loc_length * (iloc + 1); ix++)
                            for (int iy = 0; iy < dim1; iy++)
                                for (int iz = 0; iz < dim2; iz++)
                                {

                                    const int index
                                        = ix * incx + iy * incy + iz;
                                    const double alpha
                                        = (double)orbitals1.psi(color)[index]
                                          * (double)orbitals2.psi(
                                              jstate)[index];
                                    atmp[0] += alpha * cosx[ix];
                                    atmp[1] += alpha * sinx[ix];
                                    atmp[2] += alpha * cosy[iy];
                                    atmp[3] += alpha * siny[iy];
                                    atmp[4] += alpha * cosz[iz];
                                    atmp[5] += alpha * sinz[iz];
                                }

                        int ij = j * numst + i; // row i, column j

                        a[0][ij] += atmp[0];
                        a[1][ij] += atmp[1];
                        a[2][ij] += atmp[2];
                        a[3][ij] += atmp[3];
                        a[4][ij] += atmp[4];
                        a[5][ij] += atmp[5];
                    }
                }
        }
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    for (short i = 0; i < 6; i++)
    {
        mmpi.split_allreduce_sums_double(&a[i][0], n2);
        my_dscal(n2, grid.vel(), &a[i][0]);
    }

    compute_tm_.stop();
}

template <class T>
void SinCosOps<T>::computeDiag(const T& orbitals,
    VariableSizeMatrix<sparserow>& mat, const bool normalized_functions)
{
    compute_tm_.start();

    const pb::Grid& grid(orbitals.grid_);

    const int dim0 = grid.dim(0);
    const int dim1 = grid.dim(1);
    const int dim2 = grid.dim(2);

    int loc_length = dim0 / orbitals.subdivx();
    assert(loc_length > 0);
    assert(loc_length <= dim0);

    int incx = dim1 * dim2;
    int incy = dim2;

    vector<double> sinx;
    vector<double> siny;
    vector<double> sinz;
    vector<double> cosx;
    vector<double> cosy;
    vector<double> cosz;
    grid.getSinCosFunctions(sinx, siny, sinz, cosx, cosy, cosz);

    vector<vector<double>> inv_norms2;
    if (!normalized_functions)
    {
        orbitals.computeInvNorms2(inv_norms2);
    }

    // initialize sparse ordering of rows to match local overlap regions
    // This is necessary for computing correct moves in moveTo()
    mat.setupSparseRows(orbitals.getAllOverlappingGids());

    const int size = orbitals.chromatic_number();

    for (short iloc = 0; iloc < orbitals.subdivx(); iloc++)
    {
        for (short icolor = 0; icolor < size; icolor++)
        {
            int gid = orbitals.overlapping_gids_[iloc][icolor];
            if (gid != -1)
            {
                const ORBDTYPE* const psii   = orbitals.psi(icolor);
                using memory_space_type      = typename T::memory_space_type;
                unsigned int const psii_size = orbitals.getLocNumpt();
                ORBDTYPE* psii_host_view     = MemorySpace::Memory<ORBDTYPE,
                    memory_space_type>::allocate_host_view(psii_size);
                MemorySpace::Memory<ORBDTYPE, memory_space_type>::
                    copy_view_to_host(
                        const_cast<ORBDTYPE*>(psii), psii_size, psii_host_view);

                double atmp[6] = { 0., 0., 0., 0., 0., 0. };

                for (int ix = loc_length * iloc; ix < loc_length * (iloc + 1);
                     ix++)
                    for (int iy = 0; iy < dim1; iy++)
                        for (int iz = 0; iz < dim2; iz++)
                        {

                            const int index = ix * incx + iy * incy + iz;
                            const double alpha
                                = static_cast<double>(psii_host_view[index])
                                  * static_cast<double>(psii_host_view[index]);
                            atmp[0] += alpha * cosx[ix];
                            atmp[1] += alpha * sinx[ix];
                            atmp[2] += alpha * cosy[iy];
                            atmp[3] += alpha * siny[iy];
                            atmp[4] += alpha * cosz[iz];
                            atmp[5] += alpha * sinz[iz];
                        }
                MemorySpace::Memory<ORBDTYPE,
                    memory_space_type>::free_host_view(psii_host_view);
                if (!normalized_functions)
                {
                    for (int col = 0; col < 6; col++)
                        atmp[col] *= inv_norms2[iloc][icolor];
                }
                for (int col = 0; col < 6; col++)
                    mat.insertMatrixElement(gid, col, atmp[col], ADD, true);
            }
        }
    }
    /* scale data */
    mat.scale(grid.vel());
    /* gather data */
    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };

    double maxr = orbitals.getMaxR();
    DataDistribution distributor("Distributor4SinCos", maxr, myPEenv, domain);
    distributor.augmentLocalData(mat, true);

    compute_tm_.stop();
}

template class SinCosOps<LocGridOrbitals<ORBDTYPE>>;
template class SinCosOps<ExtendedGridOrbitals<ORBDTYPE>>;
