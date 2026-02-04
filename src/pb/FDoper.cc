// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "FDoper.h"
#include "FDkernels.h"
#include "memory_space.h"
#include "tools.h"

#include <iomanip>

#ifdef HAVE_MAGMA
#include "magma_v2.h"

template <typename ScalarType>
using MemoryDev = MemorySpace::Memory<ScalarType, MemorySpace::Device>;
template <typename ScalarType>
using MemoryHost = MemorySpace::Memory<ScalarType, MemorySpace::Host>;
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
}

template <class T>
void FDoper<T>::del1_2nd(
    GridFunc<T>& A, GridFunc<T>& B, const short direction) const
{
    assert(grid_.ghost_pt() > 0);

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
                u[iiz]
                    = (T)(e1 * ((double)v[iiz + incc] - (double)v[iiz - incc])
                          + e2
                                * ((double)v[iiz - incc2]
                                    - (double)v[iiz + incc2])
                          + e3
                                * ((double)v[iiz + incc3]
                                    - (double)v[iiz - incc3]));

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

                u[iiz]
                    = (T)(e1 * ((double)v[iiz + incc] - (double)v[iiz - incc])
                          + e2
                                * ((double)v[iiz - incc2]
                                    - (double)v[iiz + incc2])
                          + e3
                                * ((double)v[iiz + incc3]
                                    - (double)v[iiz - incc3])
                          + e4
                                * ((double)v[iiz - incc4]
                                    - (double)v[iiz + incc4]));

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
    if (!A.updated_boundaries()) A.trade_boundaries();

    FDkernelDel2_2nd(A.grid(), A.uu(), B.uu(), 1, MemorySpace::Host());

    B.set_updated_boundaries(0);
}

template <class T>
void FDoper<T>::del2_4th(GridFunc<T>& A, GridFunc<T>& B) const
{
    assert(grid_.ghost_pt() > 1);

    A.trade_boundaries();

    FDkernelDel2_4th(A.grid(), A.uu(), B.uu(), 1, MemorySpace::Host());

    B.set_updated_boundaries(0);
}

template <class T>
void FDoper<T>::del2_4th_withPot(
    GridFunc<T>& A, const double* const pot, T* B) const
{
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
    A.trade_boundaries();

    FDkernelDel2_6th(A.grid(), A.uu(), B.uu(), 1, MemorySpace::Host());

    B.set_updated_boundaries(0);
}

template <class T>
void FDoper<T>::del2_8th(GridFunc<T>& A, GridFunc<T>& B) const
{
    A.trade_boundaries();

    FDkernelDel2_8th(A.grid(), A.uu(), B.uu(), 1, MemorySpace::Host());

    B.set_updated_boundaries(0);
}

template <class T>
void FDoper<T>::del2_4th_Mehr(GridFunc<T>& A, GridFunc<T>& B) const
{
    A.trade_boundaries();

    FDkernelDel2_4th_Mehr(A.grid(), A.uu(0), B.uu(0), 1, MemorySpace::Host());

    B.set_updated_boundaries(0);
}

template <class T>
void FDoper<T>::smooth(GridFunc<T>& A, GridFunc<T>& B, const double alpha)
{
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
    assert(grid_.ghost_pt() > 0);

    if (!A.updated_boundaries()) A.trade_boundaries();

    rhs_4th_Mehr1_tm_.start();

    FDkernelRHS_4th_Mehr1(A.grid(), A.uu(0), B.uu(0), A.grid().ghost_pt(), 1,
        MemorySpace::Host());

    B.set_updated_boundaries(0);

    rhs_4th_Mehr1_tm_.stop();
}

template <class T>
void FDoper<T>::rhs_4th_Mehr1(GridFunc<T>& A, T* const u) const
{
    if (!A.updated_boundaries()) A.trade_boundaries();

    rhs_4th_Mehr1_tm_.start();

    FDkernelRHS_4th_Mehr1(A.grid(), A.uu(0), u, 0, 1, MemorySpace::Host());

    rhs_4th_Mehr1_tm_.stop();
}

template <class T>
void FDoper<T>::rhs_4th_Mehr2(GridFunc<T>& A, GridFunc<T>& B) const
{
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

                u0[iz] = (T)(c0 * (double)v[iiz]

                             + c1
                                   * (double)(v[iiz - incx_] + v[iiz + incx_]
                                              + v[iiz - incy_] + v[iiz + incy_]
                                              + v[iiz - 1] + v[iiz + 1])

                             + c2
                                   * (double)(v[iiz - incx_ - incy_]
                                              + v[iiz + incx_ - incy_]
                                              + v[iiz - incx_ + incy_]
                                              + v[iiz + incx_ + incy_]
                                              + v[iiz - incy_ - 1]
                                              + v[iiz - incy_ + 1]
                                              + v[iiz + incy_ - 1]
                                              + v[iiz + incy_ + 1]
                                              + v[iiz - incx_ - 1]
                                              + v[iiz - incx_ + 1]
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
