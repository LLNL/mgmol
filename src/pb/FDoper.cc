// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: FDoper.cc,v 1.22 2010/01/28 22:56:31 jeanluc Exp $
#include "FDoper.h"
#include "tools.h"
#include <iomanip>

#include "memory_space.h"

#ifdef HAVE_MAGMA
#include "magma_v2.h"

template <typename ScalarType>
using MemoryDev = MemorySpace::Memory<ScalarType, MemorySpace::Device>;
#endif

namespace pb
{

Timer FDoperInterface::del2_4th_Mehr_tm_("FDoper::del2_4th_Mehr");
Timer FDoperInterface::rhs_4th_Mehr1_tm_("FDoper::rhs_4th_Mehr1");
Timer FDoperInterface::rhs_tm_("FDoper::rhs_");
Timer FDoperInterface::del2_2nd_tm_("FDoper::del2_2nd");
Timer FDoperInterface::del2_4th_tm_("FDoper::del2_4th");
Timer FDoperInterface::del2_4th_wpot_tm_("FDoper::del2_4th_wpot");

const double inv12 = 1. / 12.;

#if 0
extern "C"
{
void lap3d_4thmehr_(const int&, const int&, const int&,
                 const int&, const int&, const int&,
                 const int&, const int&, const int&,
                 const int&, const int&, const int&,
                 const double&, const double&, const double&, const double&, 
                 const double&, const double&, const double&, 
                 const double*,
                 const double*);
}
#endif
template <class T>
FDoper<T>::FDoper(const Grid& mygrid) : grid_(mygrid)
{
    lower_order_grid_ = nullptr;
    if (grid_.active())
    {
        for (int i = 0; i < 3; i++)
        {
            assert(grid_.hgrid(i) > 1.e-8);
            inv_h_[i]  = 1. / grid_.hgrid(i);
            inv_h2_[i] = 1. / (grid_.hgrid(i) * grid_.hgrid(i));
            dim_[i]    = grid_.dim(i);
        }
        incx_ = grid_.inc(0);
        incy_ = grid_.inc(1);
        for (int i = 0; i < 3; i++)
        {
            assert(grid_.dim(i) > 0);
            assert(grid_.dim(i) < 10000);
        }
        c0mehr4_  = 16. * inv12 * (inv_h2_[0] + inv_h2_[1] + inv_h2_[2]);
        cxmehr4_  = -10. * inv12 * inv_h2_[0] + 0.125 * c0mehr4_;
        cymehr4_  = -10. * inv12 * inv_h2_[1] + 0.125 * c0mehr4_;
        czmehr4_  = -10. * inv12 * inv_h2_[2] + 0.125 * c0mehr4_;
        cxymehr4_ = -inv12 * (inv_h2_[0] + inv_h2_[1]);
        cyzmehr4_ = -inv12 * (inv_h2_[2] + inv_h2_[1]);
        cxzmehr4_ = -inv12 * (inv_h2_[0] + inv_h2_[2]);
    }
}

template <class T>
void FDoper<T>::del1_2nd(
    GridFunc<T>& A, GridFunc<T>& B, const short direction) const
{
    assert(grid_.ghost_pt() > 0);

    if (!grid_.active()) return;

    if (!A.updated_boundaries()) A.trade_boundaries();

    int iiy, iiz;
    const double e1 = 0.5 * inv_h(direction);
    const int incc  = A.grid().inc(direction);

    const T* __restrict__ v = A.uu();
    T* __restrict__ u       = B.uu();

    const int dim0 = A.grid().dim(0);
    const int dim1 = A.grid().dim(1);
    const int dim2 = A.grid().dim(2);

    const int gpt = grid_.ghost_pt();
    int iix       = gpt * incx_;

    for (int ix = 0; ix < dim0; ix++)
    {
        iiy = iix + gpt * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {
            iiz = iiy + gpt;

            for (int iz = 0; iz < dim2; iz++)
            {
                u[iiz] = e1 * (v[iiz + incc] - v[iiz - incc]);

                iiz++;
            }

            iiy += incy_;
        }

        iix += incx_;
    }

    B.set_updated_boundaries(0);
}

template <class T>
void FDoper<T>::del1_4th(
    GridFunc<T>& A, GridFunc<T>& B, const short direction) const
{
    assert(grid_.ghost_pt() > 1);
    assert(direction < 3);
    assert(direction >= 0);

    if (!grid_.active()) return;

    assert(grid_.ghost_pt() == A.grid().ghost_pt());

    if (!A.updated_boundaries()) A.trade_boundaries();

    const double e1 = (8. / 12.) * inv_h(direction);
    const double e2 = inv12 * inv_h(direction);

    const int incc  = A.grid().inc(direction);
    const int incc2 = 2 * A.grid().inc(direction);

    const int dim0 = A.grid().dim(0);
    const int dim1 = A.grid().dim(1);
    const int dim2 = A.grid().dim(2);

    const T* __restrict__ v = A.uu();
    const int gpt           = grid_.ghost_pt();
    int iix                 = gpt * incx_;

    for (int ix = 0; ix < dim0; ix++)
    {
        int iiy = iix + gpt * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {
            int iiz = iiy + gpt;

            T* __restrict__ u = B.uu() + iiy + gpt;
            assert(v != u);
            for (int iz = 0; iz < dim2; iz++)
            {
                assert(iiz < static_cast<int>(B.grid().sizeg()));
                assert((iiz + incc) < static_cast<int>(A.grid().sizeg()));
                assert((iiz + incc2) < static_cast<int>(A.grid().sizeg()));
                assert((iiz - incc) >= 0);
                assert((iiz - incc2) >= 0);

                const double v1 = (double)v[iiz + incc];
                const double v2 = (double)v[iiz - incc];
                const double v3 = (double)v[iiz + incc2];
                const double v4 = (double)v[iiz - incc2];
                u[iz]           = (T)(e1 * (v1 - v2) + e2 * (v4 - v3));

                iiz++;
            }

            iiy += incy_;
        }

        iix += incx_;
    }

    B.set_updated_boundaries(0);
}
template <class T>
void FDoper<T>::del1_6th(
    GridFunc<T>& A, GridFunc<T>& B, const short direction) const
{
    assert(direction < 3);
    assert(direction >= 0);
    assert(grid_.ghost_pt() > 2);

    if (!grid_.active()) return;

    if (!A.updated_boundaries()) A.trade_boundaries();

    const double e1 = (45. / 60.) * inv_h(direction);
    const double e2 = (9. / 60.) * inv_h(direction);
    const double e3 = (1. / 60.) * inv_h(direction);

    const int incc  = A.grid().inc(direction);
    const int incc2 = 2 * incc;
    const int incc3 = 3 * incc;

    const int dim0 = A.grid().dim(0);
    const int dim1 = A.grid().dim(1);
    const int dim2 = A.grid().dim(2);

    const T* __restrict__ v = A.uu();
    T* __restrict__ u       = B.uu();

    const int gpt = grid_.ghost_pt();
    int iix       = gpt * incx_;

    assert(v != u);

    for (int ix = 0; ix < dim0; ix++)
    {
        int iiy = iix + gpt * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {
            int iiz = iiy + gpt;

            for (int iz = 0; iz < dim2; iz++)
            {
                u[iiz] = (T)(
                    e1 * ((double)v[iiz + incc] - (double)v[iiz - incc])
                    + e2 * ((double)v[iiz - incc2] - (double)v[iiz + incc2])
                    + e3 * ((double)v[iiz + incc3] - (double)v[iiz - incc3]));

                iiz++;
            }

            iiy += incy_;
        }

        iix += incx_;
    }

    B.set_updated_boundaries(0);
}

template <class T>
void FDoper<T>::del1_8th(
    GridFunc<T>& A, GridFunc<T>& B, const short direction) const
{
    if (!grid_.active()) return;

    assert(direction < 3);
    assert(direction >= 0);
    assert(grid_.ghost_pt() > 3);

    if (!A.updated_boundaries()) A.trade_boundaries();

    const double e1 = (4. / 5.) * inv_h(direction);
    const double e2 = (1. / 5.) * inv_h(direction);
    const double e3 = (4. / 105.) * inv_h(direction);
    const double e4 = (1. / 280.) * inv_h(direction);

    const int incc  = A.grid().inc(direction);
    const int incc2 = 2 * incc;
    const int incc3 = 3 * incc;
    const int incc4 = 4 * incc;

    const int dim0 = A.grid().dim(0);
    const int dim1 = A.grid().dim(1);
    const int dim2 = A.grid().dim(2);

    const T* __restrict__ v = A.uu();
    T* __restrict__ u       = B.uu();

    const int gpt = grid_.ghost_pt();
    int iix       = gpt * incx_;

    assert(v != u);

    for (int ix = 0; ix < dim0; ix++)
    {

        int iiy = iix + gpt * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {

            int iiz = iiy + gpt;

            for (int iz = 0; iz < dim2; iz++)
            {

                u[iiz] = (T)(
                    e1 * ((double)v[iiz + incc] - (double)v[iiz - incc])
                    + e2 * ((double)v[iiz - incc2] - (double)v[iiz + incc2])
                    + e3 * ((double)v[iiz + incc3] - (double)v[iiz - incc3])
                    + e4 * ((double)v[iiz - incc4] - (double)v[iiz + incc4]));

                iiz++;
            }

            iiy += incy_;
        }

        iix += incx_;
    }

    B.set_updated_boundaries(0);
}

template <class T>
void FDoper<T>::del2_2nd(GridFunc<T>& A, GridFunc<T>& B) const
{
    if (!grid_.active()) return;

    if (!A.updated_boundaries()) A.trade_boundaries();

    del2_2nd_tm_.start();

    const T* __restrict__ v = A.uu();
    T* __restrict__ u       = B.uu();

    double cc        = inv_h2(0);
    const double c1x = -cc;

    cc               = inv_h2(1);
    const double c1y = -cc;

    cc               = inv_h2(2);
    const double c1z = -cc;
    const double c0  = -2. * (c1x + c1y + c1z);

    assert(grid_.ghost_pt() > 0);

    const int dim0 = A.grid().dim(0);
    const int dim1 = A.grid().dim(1);
    const int dim2 = A.grid().dim(2);

    const int gpt = grid_.ghost_pt();
    int iix       = gpt * incx_;

    for (int ix = 0; ix < dim0; ix++)
    {
        int iiy = iix + gpt * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {
            int iiz = iiy + gpt;

            for (int iz = 0; iz < dim2; iz++)
            {
                u[iiz] = c0 * v[iiz] + c1x * (v[iiz - incx_] + v[iiz + incx_])
                         + c1y * (v[iiz - incy_] + v[iiz + incy_])
                         + c1z * (v[iiz - 1] + v[iiz + 1]);

                iiz++;
            }

            iiy += incy_;
        }

        iix += incx_;
    }

    B.set_updated_boundaries(0);

    del2_2nd_tm_.stop();
}

template <class T>
void FDoper<T>::del2_4th(GridFunc<T>& A, GridFunc<T>& B) const
{
    if (!grid_.active()) return;

    assert(grid_.ghost_pt() > 1);

    if (!A.updated_boundaries()) A.trade_boundaries();

    del2_4th(A.grid(), A.uu(), B.uu(), 1);

    B.set_updated_boundaries(0);
}

template <class T>
void FDoper<T>::del2_4th(
    const Grid& Agrid, T* A, T* B, const size_t nfunc) const
{
    if (!grid_.active()) return;

    assert(grid_.ghost_pt() > 1);

    del2_4th_tm_.start();

    const double cc0 = inv12 * inv_h2_[0];
    const double c1x = -16. * cc0;
    const double c2x = 1. * cc0;

    const double cc1 = inv12 * inv_h2_[1];
    const double c1y = -16. * cc1;
    const double c2y = 1. * cc1;

    const double cc2 = inv12 * inv_h2_[2];
    const double c1z = -16. * cc2;
    const double c2z = 1. * cc2;

    const double c0 = -2. * (c1x + c2x + c1y + c2y + c1z + c2z);

    const int gpt = grid_.ghost_pt();

    const int incx2 = 2 * Agrid.inc(0);
    const int incy2 = 2 * Agrid.inc(1);

    int iix = gpt * incx_;

    const int dim0 = Agrid.dim(0);
    const int dim1 = Agrid.dim(1);
    const int dim2 = Agrid.dim(2);

    const size_t ng = grid_.sizeg();

#ifdef HAVE_OPENMP_OFFLOAD
    std::unique_ptr<T, void (*)(T*)> A_dev(
        MemoryDev<T>::allocate(nfunc * ng), MemoryDev<T>::free);
    MemorySpace::copy_to_dev(A, nfunc * ng, A_dev.get());

    std::unique_ptr<T, void (*)(T*)> B_dev(
        MemoryDev<T>::allocate(nfunc * ng), MemoryDev<T>::free);

    T* const A_alias = A_dev.get();
    T* B_alias       = B_dev.get();
#else
    T* const A_alias    = A;
    T* B_alias          = B;
#endif

    int incx = incx_;
    int incy = incy_;

    MGMOL_PARALLEL_FOR_COLLAPSE(4, A_alias, B_alias)
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

                    const T* __restrict__ v0   = A_alias + iiz;
                    const T* __restrict__ vmx  = A_alias + (iiz - incx);
                    const T* __restrict__ vpx  = A_alias + (iiz + incx);
                    const T* __restrict__ vmx2 = A_alias + (iiz - incx2);
                    const T* __restrict__ vpx2 = A_alias + (iiz + incx2);
                    const T* __restrict__ vmy  = A_alias + (iiz - incy);
                    const T* __restrict__ vpy  = A_alias + (iiz + incy);
                    const T* __restrict__ vmy2 = A_alias + (iiz - incy2);
                    const T* __restrict__ vpy2 = A_alias + (iiz + incy2);

                    T* __restrict__ u = B_alias + iiz;

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

#ifdef HAVE_OPENMP_OFFLOAD
    MemorySpace::copy_to_host(B_alias, ng * nfunc, B);
#endif

    del2_4th_tm_.stop();
}

template <class T>
void FDoper<T>::del2_4th_withPot(
    GridFunc<T>& A, const double* const pot, T* B) const
{
    if (!grid_.active()) return;

    assert(grid_.ghost_pt() > 1);

    if (!A.updated_boundaries()) A.trade_boundaries();

    del2_4th_wpot_tm_.start();

    const double cc0 = inv12 * inv_h2_[0];
    const double c1x = -16. * cc0;
    const double c2x = 1. * cc0;

    const double cc1 = inv12 * inv_h2_[1];
    const double c1y = -16. * cc1;
    const double c2y = 1. * cc1;

    const double cc2 = inv12 * inv_h2_[2];
    const double c1z = -16. * cc2;
    const double c2z = 1. * cc2;

    const double c0 = -2. * (c1x + c2x + c1y + c2y + c1z + c2z);

    const int gpt = grid_.ghost_pt();

    const int incx2 = 2 * A.grid().inc(0);
    const int incy2 = 2 * A.grid().inc(1);

    int iix = gpt * incx_;
    int ipx = 0;

    const int dim0 = A.dim(0);
    const int dim1 = A.dim(1);
    const int dim2 = A.dim(2);

    for (int ix = 0; ix < dim0; ix++)
    {
        int iiy = iix + gpt * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {
            int iiz = iiy + gpt;

            const T* __restrict__ v0   = A.uu(iiz);
            const T* __restrict__ vmx  = A.uu(iiz - incx_);
            const T* __restrict__ vpx  = A.uu(iiz + incx_);
            const T* __restrict__ vmx2 = A.uu(iiz - incx2);
            const T* __restrict__ vpx2 = A.uu(iiz + incx2);
            const T* __restrict__ vmy  = A.uu(iiz - incy_);
            const T* __restrict__ vpy  = A.uu(iiz + incy_);
            const T* __restrict__ vmy2 = A.uu(iiz - incy2);
            const T* __restrict__ vpy2 = A.uu(iiz + incy2);

            T* __restrict__ u = &B[ipx];

            for (int iz = 0; iz < dim2; iz++)
            {
                u[iz] = ((c0 + pot[ipx]) * (double)v0[iz]

                         + c1x * ((double)vmx[iz] + (double)vpx[iz])
                         + c1y * ((double)vmy[iz] + (double)vpy[iz])
                         + c1z * ((double)v0[iz - 1] + (double)v0[iz + 1])

                         + c2x * ((double)vmx2[iz] + (double)vpx2[iz])
                         + c2y * ((double)vmy2[iz] + (double)vpy2[iz])
                         + c2z * ((double)v0[iz - 2] + (double)v0[iz + 2]));

                ipx++;
            }

            iiy += incy_;
        }

        iix += incx_;
    }

    del2_4th_wpot_tm_.stop();
}

template <class T>
void FDoper<T>::del2_6th(GridFunc<T>& A, GridFunc<T>& B) const
{
    assert(grid_.ghost_pt() > 2);

    if (!grid_.active()) return;

    if (!A.updated_boundaries()) A.trade_boundaries();

    double cc        = (1. / 180.) * inv_h2(0);
    const double c1x = -270. * cc;
    const double c2x = 27. * cc;
    const double c3x = -2. * cc;

    cc               = (1. / 180.) * inv_h2(1);
    const double c1y = -270. * cc;
    const double c2y = 27. * cc;
    const double c3y = -2. * cc;

    cc               = (1. / 180.) * inv_h2(2);
    const double c1z = -270. * cc;
    const double c2z = 27. * cc;
    const double c3z = -2. * cc;

    const double c0
        = -2. * (c1x + c2x + c3x + c1y + c2y + c3y + c1z + c2z + c3z);

    const int incx2 = 2 * A.grid().inc(0);
    const int incy2 = 2 * A.grid().inc(1);
    const int incx3 = 3 * A.grid().inc(0);
    const int incy3 = 3 * A.grid().inc(1);

    const int dim0 = A.dim(0);
    const int dim1 = A.dim(1);
    const int dim2 = A.dim(2);

    const T* __restrict__ v = A.uu();
    T* __restrict__ u       = B.uu();
    const int gpt           = grid_.ghost_pt();

    int iix = gpt * incx_;

    for (int ix = 0; ix < dim0; ix++)
    {
        int iiy = iix + gpt * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {
            int iiz = iiy + gpt;

            for (int iz = 0; iz < dim2; iz++)
            {
                u[iiz] = (T)(
                    c0 * (double)v[iiz]

                    + c1x * ((double)v[iiz - incx_] + (double)v[iiz + incx_])
                    + c1y * ((double)v[iiz - incy_] + (double)v[iiz + incy_])
                    + c1z * ((double)v[iiz - 1] + (double)v[iiz + 1])

                    + c2x * ((double)v[iiz - incx2] + (double)v[iiz + incx2])
                    + c2y * ((double)v[iiz - incy2] + (double)v[iiz + incy2])
                    + c2z * ((double)v[iiz - 2] + (double)v[iiz + 2])

                    + c3x * ((double)v[iiz - incx3] + (double)v[iiz + incx3])
                    + c3y * ((double)v[iiz - incy3] + (double)v[iiz + incy3])
                    + c3z * ((double)v[iiz - 3] + (double)v[iiz + 3]));

                iiz++;
            }

            iiy += incy_;
        }

        iix += incx_;
    }

    B.set_updated_boundaries(0);
}
template <class T>
void FDoper<T>::del2_8th(GridFunc<T>& A, GridFunc<T>& B) const
{
    assert(grid_.ghost_pt() > 3);

    if (!grid_.active()) return;

    if (!A.updated_boundaries()) A.trade_boundaries();

    double cc        = (1. / 5040.) * inv_h2(0);
    const double c1x = -8064. * cc;
    const double c2x = 1008. * cc;
    const double c3x = -128. * cc;
    const double c4x = 9. * cc;

    cc               = (1. / 5040.) * inv_h2(1);
    const double c1y = -8064. * cc;
    const double c2y = 1008. * cc;
    const double c3y = -128. * cc;
    const double c4y = 9. * cc;

    cc               = (1. / 5040.) * inv_h2(2);
    const double c1z = -8064. * cc;
    const double c2z = 1008. * cc;
    const double c3z = -128. * cc;
    const double c4z = 9. * cc;

    const double c0 = -2.
                      * (c1x + c2x + c3x + c4x + c1y + c2y + c3y + c4y + c1z
                            + c2z + c3z + c4z);

    const int incx2 = 2 * A.grid().inc(0);
    const int incy2 = 2 * A.grid().inc(1);
    const int incx3 = 3 * A.grid().inc(0);
    const int incy3 = 3 * A.grid().inc(1);
    const int incx4 = 4 * A.grid().inc(0);
    const int incy4 = 4 * A.grid().inc(1);

    const int dim0 = A.dim(0);
    const int dim1 = A.dim(1);
    const int dim2 = A.dim(2);

    const T* __restrict__ v = A.uu();
    T* __restrict__ u       = B.uu();
    const int gpt           = grid_.ghost_pt();
    int iix                 = gpt * incx_;

    for (int ix = 0; ix < dim0; ix++)
    {
        int iiy = iix + gpt * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {
            int iiz = iiy + gpt;

            for (int iz = 0; iz < dim2; iz++)
            {
                u[iiz] = (T)(
                    c0 * (double)v[iiz]

                    + c1x * ((double)v[iiz - incx_] + (double)v[iiz + incx_])
                    + c1y * ((double)v[iiz - incy_] + (double)v[iiz + incy_])
                    + c1z * ((double)v[iiz - 1] + (double)v[iiz + 1])

                    + c2x * ((double)v[iiz - incx2] + (double)v[iiz + incx2])
                    + c2y * ((double)v[iiz - incy2] + (double)v[iiz + incy2])
                    + c2z * ((double)v[iiz - 2] + (double)v[iiz + 2])

                    + c3x * ((double)v[iiz - incx3] + (double)v[iiz + incx3])
                    + c3y * ((double)v[iiz - incy3] + (double)v[iiz + incy3])
                    + c3z * ((double)v[iiz - 3] + (double)v[iiz + 3])

                    + c4x * ((double)v[iiz - incx4] + (double)v[iiz + incx4])
                    + c4y * ((double)v[iiz - incy4] + (double)v[iiz + incy4])
                    + c4z * ((double)v[iiz - 4] + (double)v[iiz + 4]));

                iiz++;
            }

            iiy += incy_;
        }

        iix += incx_;
    }

    B.set_updated_boundaries(0);
}
template <class T>
void FDoper<T>::del2_4th_Mehr(GridFunc<T>& A, GridFunc<T>& B) const
{
    if (!grid_.active()) return;

    if (!A.updated_boundaries()) A.trade_boundaries();

    del2_4th_Mehr_tm_.start();

    assert(inv_h2(0) > 1.e-12);
    assert(inv_h2(1) > 1.e-12);
    assert(inv_h2(2) > 1.e-12);
    assert(ghosts() > 0);

    const int shift = grid_.ghost_pt();
    const int dim0  = grid_.dim(0);
    const int dim1  = grid_.dim(1);
    const int dim2  = grid_.dim(2);

#if 0
    int ifirst0=1+shift;
    int ilast0=dim0+shift;
    int ifirst1=1+shift;
    int ilast1=dim1+shift;
    int ifirst2=1+shift;
    int ilast2=dim2+shift;
    int ilo0=1;
    int ihi0=dim0+2*shift;
    int ilo1=1;
    int ihi1=dim1+2*shift;
    int ilo2=1;
    int ihi2=dim2+2*shift;
    lap3d_4thmehr_(ifirst2,ilast2,ifirst1,ilast1,ifirst0,ilast0,
                   ilo2,ihi2,ilo1,ihi1,ilo0,ihi0,
                   c0mehr4_,czmehr4_,cymehr4_,cxmehr4_,cyzmehr4_,cxymehr4_,cxzmehr4_,
                   A.uu(0),B.uu(0));
#else
    const int iix0      = shift * incx_;

    const T* const v = A.uu(0);
    T* u             = B.uu(0);

    for (int ix = 0; ix < dim0; ix++)
    {

        int iix = iix0 + ix * incx_;
        int iiy = iix + shift * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {

            const int iiz = iiy + shift;

            T* const u0          = u + iiz;
            const T* const v0    = v + iiz;
            const T* const vmx   = v0 - incx_;
            const T* const vpx   = v0 + incx_;
            const T* const vmxmy = vmx - incy_;
            const T* const vpxmy = vpx - incy_;
            const T* const vmy   = v0 - incy_;
            const T* const vpy   = v0 + incy_;
            const T* const vmxpy = vmx + incy_;
            const T* const vpxpy = vpx + incy_;

            for (int iz = 0; iz < dim2; iz++)
            {

                u0[iz] = (T)(c0mehr4_ * (double)v0[iz]

                             + czmehr4_ * (double)(v0[iz - 1] + v0[iz + 1])
                             + cymehr4_ * (double)(vmy[iz] + vpy[iz])
                             + cxmehr4_ * (double)(vmx[iz] + vpx[iz])

                             + cxzmehr4_
                                   * (double)(vmx[iz - 1] + vmx[iz + 1]
                                              + vpx[iz - 1] + vpx[iz + 1])
                             + cyzmehr4_
                                   * (double)(vmy[iz - 1] + vmy[iz + 1]
                                              + vpy[iz - 1] + vpy[iz + 1])
                             + cxymehr4_
                                   * (double)(vmxmy[iz] + vpxmy[iz] + vmxpy[iz]
                                              + vpxpy[iz]));
            }

            iiy += incy_;
        }
    }
#endif
    B.set_updated_boundaries(0);

    del2_4th_Mehr_tm_.stop();
}

template <class T>
void FDoper<T>::smooth(GridFunc<T>& A, GridFunc<T>& B, const double alpha)
{
    if (!grid_.active()) return;

    if (!A.updated_boundaries()) A.trade_boundaries();

    const double c0 = alpha;
    const double c1 = (1. - alpha) / 6.;

    assert(grid_.ghost_pt() > 0);

    T* const u       = B.uu();
    const T* const v = A.uu();
    const int shift  = grid_.ghost_pt();
    const int iix0   = shift * incx_;

    const int dim0 = dim(0);
    const int dim1 = dim(1);
    const int dim2 = dim(2);
    for (int ix = 0; ix < dim0; ix++)
    {

        int iix = iix0 + ix * incx_;
        int iiy = iix + shift * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {

            int iiz = iiy + shift;

            for (int iz = 0; iz < dim2; iz++)
            {

                u[iiz] = (T)(c0 * (double)v[iiz]

                             + c1
                                   * (double)(v[iiz - incx_] + v[iiz + incx_]
                                              + v[iiz - incy_] + v[iiz + incy_]
                                              + v[iiz - 1] + v[iiz + 1]));

                iiz++;
            }

            iiy += incy_;
        }
    }

    B.set_updated_boundaries(0);
}

// void rhs(GridFunc<T> &A, GridFunc<T> &B)const{ B=A; };
template <class T>
void FDoper<T>::rhs_4th_Mehr1(GridFunc<T>& A, GridFunc<T>& B) const
{
    if (!grid_.active()) return;

    assert(grid_.ghost_pt() > 0);

    if (!A.updated_boundaries()) A.trade_boundaries();

    rhs_4th_Mehr1_tm_.start();

    const double c0 = 0.5;
    const double c1 = inv12;

    const int shift = grid_.ghost_pt();
    const int iix0  = shift * incx_;

    const T* const v = A.uu(0);
    T* u             = B.uu(0);

    const int dim0 = dim(0);
    const int dim1 = dim(1);
    const int dim2 = dim(2);
    for (int ix = 0; ix < dim0; ix++)
    {

        int iix = iix0 + ix * incx_;
        int iiy = iix + shift * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {

            int iiz = iiy + shift;

            T* const u0        = u + iiz;
            const T* const v0  = v + iiz;
            const T* const vmx = v0 - incx_;
            const T* const vpx = v0 + incx_;
            const T* const vmy = v0 - incy_;
            const T* const vpy = v0 + incy_;

            for (int iz = 0; iz < dim2; iz++)
            {

                u0[iz]
                    = (T)(c0 * (double)v0[iz]

                          + c1
                                * (double)(vmx[iz] + vpx[iz] + vmy[iz] + vpy[iz]
                                           + v0[iz - 1] + v0[iz + 1]));
            }

            iiy += incy_;
        }
    }

    B.set_updated_boundaries(0);

    rhs_4th_Mehr1_tm_.stop();
}

template <class T>
void FDoper<T>::rhs_4th_Mehr1(GridFunc<T>& A, T* const u) const
{
    if (!grid_.active()) return;

    if (!A.updated_boundaries()) A.trade_boundaries();

    rhs_4th_Mehr1_tm_.start();

    const double c0 = 0.5;
    const double c1 = inv12;

    assert(grid_.ghost_pt() > 0);

    const int shift = grid_.ghost_pt();
    const int iix0  = shift * incx_;

    const T* const v = A.uu(0);

    const int dim0  = dim(0);
    const int dim1  = dim(1);
    const int dim2  = dim(2);
    const int dim12 = dim1 * dim2;
    for (int ix = 0; ix < dim0; ix++)
    {

        int iix = iix0 + ix * incx_;
        int iiy = iix + shift * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {

            int iiz = iiy + shift;

            T* const u0        = u + ix * dim12 + iy * dim2;
            const T* const v0  = v + iiz;
            const T* const vmx = v0 - incx_;
            const T* const vpx = v0 + incx_;
            const T* const vmy = v0 - incy_;
            const T* const vpy = v0 + incy_;

            for (int iz = 0; iz < dim2; iz++)
            {

                u0[iz]
                    = (T)(c0 * (double)v0[iz]

                          + c1
                                * (double)(vmx[iz] + vpx[iz] + vmy[iz] + vpy[iz]
                                           + v0[iz - 1] + v0[iz + 1]));
            }

            iiy += incy_;
        }
    }

    rhs_4th_Mehr1_tm_.stop();
}

template <class T>
void FDoper<T>::rhs_4th_Mehr2(GridFunc<T>& A, GridFunc<T>& B) const
{
    if (!grid_.active()) return;

    int iiy, iiz;

    const double c0 = 2. / 3.;
    const double c1 = 1. / 36.;
    const double c2 = 1. / 72.;

    assert(grid_.ghost_pt() > 1);

    if (!A.updated_boundaries()) A.trade_boundaries();

    const T* __restrict__ v = A.uu();
    T* __restrict__ u       = B.uu();

    const int gpt = grid_.ghost_pt();
    int iix       = gpt * incx_;

    const int dim0 = dim(0);
    const int dim1 = dim(1);
    const int dim2 = dim(2);

    for (int ix = 0; ix < dim0; ix++)
    {

        iiy = iix + gpt * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {

            iiz = iiy + gpt;

            for (int iz = 0; iz < dim2; iz++)
            {

                u[iiz]
                    = (T)c0 * (double)v[iiz]

                      + c1
                            * (double)(v[iiz - incx_] + v[iiz + incx_]
                                       + v[iiz - incy_] + v[iiz + incy_]
                                       + v[iiz - 1] + v[iiz + 1])

                      + c2
                            * (double)(v[iiz - incx_ - incy_]
                                       + v[iiz + incx_ - incy_]
                                       + v[iiz - incx_ + incy_]
                                       + v[iiz + incx_ + incy_]
                                       + v[iiz - incy_ - 1] + v[iiz - incy_ + 1]
                                       + v[iiz + incy_ - 1] + v[iiz + incy_ + 1]
                                       + v[iiz - incx_ - 1] + v[iiz - incx_ + 1]
                                       + v[iiz + incx_ - 1]
                                       + v[iiz + incx_ + 1]);
                iiz++;
            }

            iiy += incy_;
        }

        iix += incx_;
    }

    B.set_updated_boundaries(0);
}
template <class T>
void FDoper<T>::rhs_4th_Mehr2(GridFunc<T>& A, T* const u) const
{
    if (!grid_.active()) return;

    int iiy, iiz;

    const double c0 = 2. / 3.;
    const double c1 = 1. / 36.;
    const double c2 = 1. / 72.;

    assert(grid_.ghost_pt() > 1);

    if (!A.updated_boundaries()) A.trade_boundaries();

    const T* __restrict__ v = A.uu();

    const int gpt = grid_.ghost_pt();
    int iix       = gpt * incx_;

    const int dim0  = dim(0);
    const int dim1  = dim(1);
    const int dim2  = dim(2);
    const int dim12 = dim1 * dim2;

    for (int ix = 0; ix < dim0; ix++)
    {

        iiy = iix + gpt * incy_;

        for (int iy = 0; iy < dim1; iy++)
        {

            iiz = iiy + gpt;

            T* __restrict__ u0 = u + ix * dim12 + iy * dim2;
            for (int iz = 0; iz < dim2; iz++)
            {

                u0[iz] = (T)(
                    c0 * (double)v[iiz]

                    + c1
                          * (double)(v[iiz - incx_] + v[iiz + incx_]
                                     + v[iiz - incy_] + v[iiz + incy_]
                                     + v[iiz - 1] + v[iiz + 1])

                    + c2
                          * (double)(v[iiz - incx_ - incy_]
                                     + v[iiz + incx_ - incy_]
                                     + v[iiz - incx_ + incy_]
                                     + v[iiz + incx_ + incy_]
                                     + v[iiz - incy_ - 1] + v[iiz - incy_ + 1]
                                     + v[iiz + incy_ - 1] + v[iiz + incy_ + 1]
                                     + v[iiz - incx_ - 1] + v[iiz - incx_ + 1]
                                     + v[iiz + incx_ - 1]
                                     + v[iiz + incx_ + 1]));
                iiz++;
            }

            iiy += incy_;
        }

        iix += incx_;
    }
}
template class FDoper<double>;
template class FDoper<float>;
} // namespace pb
