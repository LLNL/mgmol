// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "XYZOps.h"

#include "ExtendedGridOrbitals.h"
#include "FunctionsPacking.h"
#include "LocGridOrbitals.h"
#include "MGmol_MPI.h"
#include "memory_space.h"

#ifdef HAVE_MAGMA
using memory_space_type = MemorySpace::Device;
#else
using memory_space_type = MemorySpace::Host;
#endif

// compute 3 symmetric matrices with expectations of X, Y and Z operators
template <class OrbitalsType>
void XYZOps<OrbitalsType>::compute(
    const OrbitalsType& orbitals, std::vector<std::vector<double>>& a)
{
    assert(a.size() == 3);

    compute_tm_.start();

    auto first_value       = orbitals.getPsi(0)[0];
    using OrbitalsDataType = decltype(first_value);

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

    std::vector<double> x, y, z;
    grid.getXYZFunctions(x, y, z);

    const int msize = orbitals.chromatic_number();

    const int ld                = orbitals.getLda();
    unsigned int const size_psi = msize * ld;
    OrbitalsDataType* psi_view  = MemorySpace::Memory<OrbitalsDataType,
        memory_space_type>::allocate_host_view(size_psi);
    MemorySpace::Memory<OrbitalsDataType, memory_space_type>::copy_view_to_host(
        orbitals.psi(0), size_psi, psi_view);

    for (short iloc = 0; iloc < orbitals.subdivx(); iloc++)
    {
        for (int icolor = 0; icolor < msize; icolor++)
        {
            const int i = orbitals.overlapping_gids_[iloc][icolor];
            if (i != -1)
            {
                const OrbitalsDataType* const ppsii = psi_view + ld * icolor;
                for (int jstate = 0; jstate <= icolor; jstate++)
                {
                    const int j = orbitals.overlapping_gids_[iloc][jstate];
                    if (j != -1)
                    {
                        const OrbitalsDataType* const ppsij
                            = psi_view + ld * jstate;

                        double atmp[6]  = { 0., 0., 0., 0., 0., 0. };
                        const int ixend = loc_length * (iloc + 1);

                        // loop over patch
                        for (int ix = loc_length * iloc; ix < ixend; ix++)
                        {
                            const double vx   = x[ix];
                            const int offsetx = ix * incx;

                            for (int iy = 0; iy < dim1; iy++)
                            {
                                const double vy  = y[iy];
                                const int offset = offsetx + iy * incy;

                                for (int iz = 0; iz < dim2; iz++)
                                {
                                    const int index    = offset + iz;
                                    const double alpha = (double)ppsij[index]
                                                         * (double)ppsii[index];
                                    atmp[0] += alpha * vx;
                                    atmp[1] += alpha * vy;
                                    atmp[2] += alpha * z[iz];
                                }
                            }
                        }
                        const int ji = j * numst + i;
                        const int ij = i * numst + j;
                        a[0][ji]     = a[0][ij] += atmp[0];
                        a[1][ji]     = a[1][ij] += atmp[1];
                        a[2][ji]     = a[2][ij] += atmp[2];
                    }
                }
            }
        }
    }

    MemorySpace::Memory<OrbitalsDataType, memory_space_type>::free_host_view(
        psi_view);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    for (short i = 0; i < 6; i++)
    {
        mmpi.split_allreduce_sums_double(&a[i][0], n2);
        my_dscal(n2, grid.vel(), &a[i][0]);
    }

    compute_tm_.stop();
}

template class XYZOps<LocGridOrbitals<ORBDTYPE>>;
template class XYZOps<ExtendedGridOrbitals<ORBDTYPE>>;
