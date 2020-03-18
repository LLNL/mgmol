// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "GridMaskMult.h"

template <typename T>
void GridMaskMult::applyPrivate(
    T* u, const unsigned short level, const unsigned short iloc)
{
    // early return if nothing to do
    if (GridMask::maskIs1(level, iloc)) return;

    if (GridMask::maskIs0(level, iloc))
    {
        GridMask::setZero(u, level, iloc);
    }
    else if (GridMask::maskIs2(level, iloc))
    {
        GridMask::multiplyByMask(u, level, iloc);
    }
}

template <typename T>
void GridMaskMult::applyPrivate(
    pb::GridFunc<T>& gu, const unsigned short level, const unsigned short iloc)
{
    gu.set_updated_boundaries(false);

    // early return if nothing to do
    if (GridMask::maskIs1(level, iloc)) return;

    const int subdimx = GridMask::subdim0(level);

    const int shift = gu.ghost_pt();
    const int incx1 = gu.grid().inc(0);

    const int offset = (shift + subdimx * iloc) * incx1;
    T* const pu      = gu.uu() + offset;

    if (GridMask::maskIs0(level, iloc))
    {

        GridMask::setZero(gu, level, iloc);
    }
    else if (GridMask::maskIs2(level, iloc))
    {

        const int incy1 = gu.grid().inc(1);

        const int dim1  = gu.grid().dim(1);
        const int dim2  = gu.grid().dim(2);
        const int incx2 = dim2 * dim1;
        const int incy2 = dim2;

        const lmasktype* const maskptr = GridMask::lmask(level, iloc);

        int ix1 = shift + shift * incy1;
        int ix2 = 0;
        for (int ix = 0; ix < subdimx; ix++)
        {
            int iy1 = ix1;
            int iy2 = ix2;
            for (int iy = 0; iy < dim1; iy++)
            {
                T* const ppu                 = pu + iy1;
                const lmasktype* const pmask = maskptr + iy2;
                for (int iz = 0; iz < dim2; iz++)
                {
                    // assert( iy2+iz< gu.grid().size()/subdivx_ );
                    assert(offset + iy1 + iz
                           < static_cast<int>(gu.grid().sizeg()));
                    ppu[iz] *= (T)pmask[iz];
                }
                iy1 += incy1;
                iy2 += incy2;
            }
            ix1 += incx1;
            ix2 += incx2;
        }
    }
}

// explicit instantiations
template void GridMaskMult::applyPrivate(pb::GridFunc<double>& gu,
    const unsigned short level, const unsigned short iloc);
template void GridMaskMult::applyPrivate(pb::GridFunc<float>& gu,
    const unsigned short level, const unsigned short iloc);
template void GridMaskMult::applyPrivate(
    float* u, const unsigned short level, const unsigned short iloc);
template void GridMaskMult::applyPrivate(
    double* u, const unsigned short level, const unsigned short iloc);
