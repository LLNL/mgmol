// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "GridFuncVector.h"
#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "mputils.h"

#include <cassert>
#include <set>

namespace pb
{
Timer GridFuncVectorInterface::trade_bc_tm_("GridFuncVector::trade_bc");
Timer GridFuncVectorInterface::trade_bc_colors_tm_(
    "GridFuncVector::trade_bc_colors");
Timer GridFuncVectorInterface::prod_tm_("GridFuncVector::prod");
Timer GridFuncVectorInterface::finishExchangeNorthSouth_tm_(
    "GridFuncVector::finishExNorthSouth");
Timer GridFuncVectorInterface::finishExchangeUpDown_tm_(
    "GridFuncVector::finishExUpDown");
Timer GridFuncVectorInterface::finishExchangeEastWest_tm_(
    "GridFuncVector::finishExEastWest");
Timer GridFuncVectorInterface::wait_north_south_tm_("GridFuncVector::waitNS");
Timer GridFuncVectorInterface::wait_up_down_tm_("GridFuncVector::waitUD");
Timer GridFuncVectorInterface::wait_east_west_tm_("GridFuncVector::waitEW");

template <typename ScalarType>
std::vector<std::vector<ScalarType>> GridFuncVector<ScalarType>::comm_buf1_;
template <typename ScalarType>
std::vector<std::vector<ScalarType>> GridFuncVector<ScalarType>::comm_buf2_;
template <typename ScalarType>
std::vector<std::vector<ScalarType>> GridFuncVector<ScalarType>::comm_buf3_;
template <typename ScalarType>
std::vector<std::vector<ScalarType>> GridFuncVector<ScalarType>::comm_buf4_;

template <typename ScalarType>
void GridFuncVector<ScalarType>::allocate(const int n)
{
    functions_.resize(n);

    // allocate memory for all GridFunc
    int alloc_size = grid_.sizeg();
    memory_.reset(new ScalarType[n * alloc_size]);

    memset(memory_.get(), 0, n * alloc_size * sizeof(ScalarType));

    for (int i = 0; i < n; i++)
    {
        ScalarType* alloc = memory_.get() + i * alloc_size;
        functions_[i]
            = new GridFunc<ScalarType>(grid_, bc_[0], bc_[1], bc_[2], alloc);
    }

    // jlf, 8/6/2020: one should be able to set this flag to true
    // but we may need to fix a few things for that to work
    updated_boundaries_ = false;
}

template <typename ScalarType>
void GridFuncVector<ScalarType>::setup()
{
    const int mytask = grid_.mype_env().mytask();
    north_ = (bc_[1] == 1 || (grid_.mype_env().mpi_neighbors(NORTH) > mytask));
    south_ = (bc_[1] == 1 || (grid_.mype_env().mpi_neighbors(SOUTH) < mytask));
    up_    = (bc_[2] == 1 || (grid_.mype_env().mpi_neighbors(UP) > mytask));
    down_  = (bc_[2] == 1 || (grid_.mype_env().mpi_neighbors(DOWN) < mytask));
    east_  = (bc_[0] == 1 || (grid_.mype_env().mpi_neighbors(EAST) > mytask));
    west_  = (bc_[0] == 1 || (grid_.mype_env().mpi_neighbors(WEST) < mytask));

    nfunc_ = (int)gid_[0].size();
    MPI_Allreduce(&nfunc_, &nfunc_max_global_, 1, MPI_INT, MPI_MAX, comm_);
    nsubdivx_ = (short)gid_.size();

    remote_gids_.resize(6);

    for (short iloc = 0; iloc < nsubdivx_; iloc++)
        for (short k = 0; k < nfunc_; k++)
            if (gid_[iloc][k] >= 0)
                gid2lid_.insert(std::pair<int, short>(gid_[iloc][k], k));

    nghosts_ = grid_.ghost_pt();

    dimx_ = grid_.dim(0);
    dimy_ = grid_.dim(1);
    dimz_ = grid_.dim(2);

    incx_ = grid_.inc(0);
    incy_ = grid_.inc(1);

    dimxy_ = (dimy_ + 2 * nghosts_) * dimx_;
    incxy_ = dimxy_ / nsubdivx_;

    east_west_size_   = nghosts_ * incx_;
    south_north_size_ = nghosts_ * dimz_ * dimx_ + nsubdivx_;
    up_down_size_     = nghosts_ * dimxy_ + nsubdivx_;

    assert(dimx_ >= nghosts_);
    assert(dimy_ >= nghosts_);
    assert(dimz_ >= nghosts_);
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::prod(
    GridFuncVector<ScalarType>& A, const GridFunc<double>& B)
{
    assert(A.grid_.sizeg() == grid_.sizeg());
    assert(B.grid().sizeg() == grid_.sizeg());
    assert(A.size() == size());

    prod_tm_.start();

    const int bsize = 64;
    const int ng    = grid_.sizeg();
    const int nb    = ng / bsize;
    const int nf    = (int)A.size();

    // loop over blocks
    for (int ib = 0; ib <= nb; ib++)
    {

        const int ibstart             = ib * bsize;
        const double* __restrict__ v2 = B.uu(ibstart);
        const int npt                 = (ib < nb) ? bsize : ng - ibstart;
        assert(npt >= 0);

        if (npt > 0)
            for (int j = 0; j < nf; j++)
            {

                ScalarType* __restrict__ pu = functions_[j]->uu(ibstart);
                const ScalarType* __restrict__ v1
                    = A.functions_[j]->uu(ibstart);
                for (int i = 0; i < npt; i++)
                {
                    pu[i] = (ScalarType)(v1[i] * v2[i]);
                }
            }
    }
    for (int j = 0; j < nf; j++)
    {
        functions_[j]->set_updated_boundaries(
            A.functions_[j]->updated_boundaries() && B.updated_boundaries());
    }

    prod_tm_.stop();
}

template <typename ScalarType>
void GridFuncVector<ScalarType>::prod(
    GridFuncVector<ScalarType>& A, const GridFunc<float>& B)
{
    assert(A.grid_.sizeg() == grid_.sizeg());
    assert(B.grid().sizeg() == grid_.sizeg());
    assert(A.size() == size());

    prod_tm_.start();

    const int bsize = 64;
    const int ng    = grid_.sizeg();
    const int nb    = ng / bsize;
    const int nf    = (int)A.size();

    // loop over blocks
    for (int ib = 0; ib <= nb; ib++)
    {

        const int ibstart            = ib * bsize;
        const float* __restrict__ v2 = B.uu(ibstart);
        const int npt                = (ib < nb) ? bsize : ng - ibstart;
        assert(npt >= 0);

        if (npt > 0)
            for (int j = 0; j < nf; j++)
            {

                ScalarType* __restrict__ pu = functions_[j]->uu(ibstart);
                const ScalarType* __restrict__ v1
                    = A.functions_[j]->uu(ibstart);
                for (int i = 0; i < npt; i++)
                {
                    pu[i] = (ScalarType)(v1[i] * v2[i]);
                }
            }
    }
    for (int j = 0; j < nf; j++)
    {
        functions_[j]->set_updated_boundaries(
            A.functions_[j]->updated_boundaries() && B.updated_boundaries());
    }

    prod_tm_.stop();
}

template <typename ScalarType>
void GridFuncVector<ScalarType>::wait_north_south()
{
    if (grid_.mype_env().n_mpi_task(1) == 1) return;

    wait_north_south_tm_.start();

    if (south_)
    {
        MPI_Wait(req_north_south_ + 2, MPI_STATUS_IGNORE);
        MPI_Wait(req_north_south_ + 1, MPI_STATUS_IGNORE);
    }
    if (north_)
    {
        MPI_Wait(req_north_south_, MPI_STATUS_IGNORE);
        MPI_Wait(req_north_south_ + 3, MPI_STATUS_IGNORE);
    }

    wait_north_south_tm_.stop();
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::wait_east_west()
{
    if (grid_.mype_env().n_mpi_task(0) == 1) return;

    wait_east_west_tm_.start();

    if (west_)
    {
        MPI_Wait(req_east_west_ + 1, MPI_STATUS_IGNORE);
        MPI_Wait(req_east_west_ + 2, MPI_STATUS_IGNORE);
    }
    if (east_)
    {
        MPI_Wait(req_east_west_, MPI_STATUS_IGNORE);
        MPI_Wait(req_east_west_ + 3, MPI_STATUS_IGNORE);
    }

    wait_east_west_tm_.stop();
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::wait_up_down()
{
    if (grid_.mype_env().n_mpi_task(2) == 1) return;

    wait_up_down_tm_.start();

    if (up_)
    {
        MPI_Wait(req_up_down_, MPI_STATUS_IGNORE);
        MPI_Wait(req_up_down_ + 3, MPI_STATUS_IGNORE);
    }
    if (down_)
    {
        MPI_Wait(req_up_down_ + 2, MPI_STATUS_IGNORE);
        MPI_Wait(req_up_down_ + 1, MPI_STATUS_IGNORE);
    }

    wait_up_down_tm_.stop();
}

template <typename ScalarType>
void GridFuncVector<ScalarType>::allocate_buffers(const int nfunc)
{
    static int last_nfunc = 0;

    if (nfunc <= last_nfunc) return;

    if (grid_.mype_env().n_mpi_tasks() > 1)
    {
        int size_max = nfunc * nghosts_ * incx_;
        size_max     = std::max(size_max, nfunc * nghosts_ * dimxy_);
        size_max     = std::max(size_max, nfunc * nghosts_ * dimx_ * incy_);
        size_max += nfunc * nsubdivx_; // to pack gids
        size_max += 1; // to pack number of functions (data) in buffer

        comm_buf1_.resize(3);
        comm_buf2_.resize(3);
        comm_buf3_.resize(3);
        comm_buf4_.resize(3);

        for (short i = 0; i < 3; i++)
        {
            comm_buf1_[i].resize(size_max);
            comm_buf2_[i].resize(size_max);
            comm_buf3_[i].resize(size_max);
            comm_buf4_[i].resize(size_max);
        }

        last_nfunc = nfunc;
    }
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::initiateNorthSouthComm(
    const int begin_color, const int end_color)
{
    std::vector<ScalarType>& comm_buf1(comm_buf1_[1]);
    std::vector<ScalarType>& comm_buf2(comm_buf2_[1]);
    std::vector<ScalarType>& comm_buf3(comm_buf3_[1]);
    std::vector<ScalarType>& comm_buf4(comm_buf4_[1]);

    const int ncolors = end_color - begin_color;

    const size_t sizebuffer = 1 + nfunc_max_global_ * south_north_size_;
    const size_t sdimz      = dimz_ * sizeof(ScalarType);
    if (north_)
        grid_.mype_env().Irecv(
            &comm_buf4[0], sizebuffer, NORTH, &req_north_south_[3]);
    if (south_)
        grid_.mype_env().Irecv(
            &comm_buf3[0], sizebuffer, SOUTH, &req_north_south_[1]);

    const int imin = nghosts_ * incx_;
    const int iinc = incx_ * dimx_ / nsubdivx_;

    const int jmax = nghosts_ * incy_;
    const int ymax = dimy_ * grid_.inc(1);

    if (south_)
    {
        ScalarType* buf2_ptr = &comm_buf2[0];
        // first element will tell how many functions (data) are in buffer
        *buf2_ptr = (ScalarType)ncolors;
        buf2_ptr++;
        // pack data
        for (int color = begin_color; color < end_color; color++)
        {
            for (short iloc = 0; iloc < nsubdivx_; iloc++)
            {
                const ScalarType* uus
                    = functions_[color]->uu(nghosts_ * (incy_ + 1));
                *buf2_ptr = (ScalarType)gid_[iloc][color];
                buf2_ptr++;
                for (int j = 0; j < jmax; j += incy_)
                    for (int i = imin + iloc * iinc;
                         i < imin + (iloc + 1) * iinc; i += incx_)
                    {
                        memcpy(buf2_ptr, &uus[i + j], sdimz);
                        buf2_ptr += dimz_;
                    }
            }
        }
        grid_.mype_env().Isend(&comm_buf2[0], 1 + ncolors * south_north_size_,
            SOUTH, &req_north_south_[2]);
    }

    if (north_)
    {
        ScalarType* buf1_ptr = &comm_buf1[0];
        // first element will tell how many functions (data) are in buffer
        *buf1_ptr = (ScalarType)ncolors;
        buf1_ptr++;
        // pack data
        for (int color = begin_color; color < end_color; color++)
        {
            for (short iloc = 0; iloc < nsubdivx_; iloc++)
            {
                const ScalarType* uus = functions_[color]->uu(nghosts_ + ymax);
                *buf1_ptr             = (ScalarType)gid_[iloc][color];
                buf1_ptr++;
                for (int j = 0; j < jmax; j += incy_)
                    for (int i = imin + iloc * iinc;
                         i < imin + (iloc + 1) * iinc; i += incx_)
                    {
                        memcpy(buf1_ptr, &uus[i + j], sdimz);
                        buf1_ptr += dimz_;
                    }
            }
        }
        grid_.mype_env().Isend(&comm_buf1[0], 1 + ncolors * south_north_size_,
            NORTH, &req_north_south_[0]);
    }
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::finishNorthSouthComm()
{
    finishExchangeNorthSouth_tm_.start();

    std::vector<ScalarType>& comm_buf3(comm_buf3_[1]);
    std::vector<ScalarType>& comm_buf4(comm_buf4_[1]);

    const int ymax     = dimy_ * grid_.inc(1);
    const size_t sdimz = dimz_ * sizeof(ScalarType);

    if (grid_.mype_env().n_mpi_task(1) > 1)
    {

        const int imin = nghosts_ * incx_;
        const int iinc = incx_ * dimx_ / nsubdivx_;
        const int jmax = nghosts_ * incy_;

        if (north_)
        {
            ScalarType* buf4_ptr   = &comm_buf4[0];
            const int nremote_func = (int)(*buf4_ptr);
            buf4_ptr++;

            for (int k = 0; k < nremote_func; k++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const int gid = (int)(*buf4_ptr);
                    buf4_ptr++;

                    if (gid >= 0)
                    {
                        std::map<int, short>::const_iterator ilid
                            = gid2lid_.find(gid);
                        if (ilid != gid2lid_.end())
                        {
                            const short lid = ilid->second;
                            if (gid == gid_[iloc][lid])
                            {
                                ScalarType* uus = functions_[lid]->uu(
                                    nghosts_ * (incy_ + 1) + ymax);
                                for (int j = 0; j < jmax; j += incy_)
                                    for (int i = imin + iloc * iinc;
                                         i < imin + (iloc + 1) * iinc;
                                         i += incx_)
                                    {
                                        memcpy(&uus[i + j], buf4_ptr, sdimz);
                                        buf4_ptr += dimz_;
                                    }
                            }
                            else
                            {
                                for (int j = 0; j < jmax; j += incy_)
                                    for (int i = imin + iloc * iinc;
                                         i < imin + (iloc + 1) * iinc;
                                         i += incx_)
                                        buf4_ptr += dimz_;
                            }
                        }
                        else
                        {
                            for (int j = 0; j < jmax; j += incy_)
                                for (int i = imin + iloc * iinc;
                                     i < imin + (iloc + 1) * iinc; i += incx_)
                                    buf4_ptr += dimz_;
                        }
                    }
                    else
                    {
                        for (int j = 0; j < jmax; j += incy_)
                            for (int i = imin + iloc * iinc;
                                 i < imin + (iloc + 1) * iinc; i += incx_)
                            {
                                buf4_ptr += dimz_;
                            }
                    }
                }
            }
            assert((buf4_ptr - &comm_buf4[0])
                   <= static_cast<int>(comm_buf4.size()));
        }

        if (south_)
        {
            ScalarType* buf3_ptr   = &comm_buf3[0];
            const int nremote_func = (int)(*buf3_ptr);
            buf3_ptr++;

            for (int k = 0; k < nremote_func; k++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const int gid = (int)(*buf3_ptr);
                    buf3_ptr++;
                    if (gid >= 0)
                    {
                        std::map<int, short>::const_iterator ilid
                            = gid2lid_.find(gid);
                        if (ilid != gid2lid_.end())
                        {
                            const short lid = ilid->second;
                            if (gid == gid_[iloc][lid])
                            {
                                ScalarType* uus = functions_[lid]->uu(nghosts_);
                                for (int j = 0; j < jmax; j += incy_)
                                    for (int i = imin + iloc * iinc;
                                         i < imin + (iloc + 1) * iinc;
                                         i += incx_)
                                    {
                                        memcpy(&uus[i + j], buf3_ptr, sdimz);
                                        buf3_ptr += dimz_;
                                    }
                            }
                            else
                            {
                                for (int j = 0; j < jmax; j += incy_)
                                    for (int i = imin + iloc * iinc;
                                         i < imin + (iloc + 1) * iinc;
                                         i += incx_)
                                        buf3_ptr += dimz_;
                            }
                        }
                        else
                        {
                            for (int j = 0; j < jmax; j += incy_)
                                for (int i = imin + iloc * iinc;
                                     i < imin + (iloc + 1) * iinc; i += incx_)
                                    buf3_ptr += dimz_;
                        }
                    }
                    else
                    {
                        for (int j = 0; j < jmax; j += incy_)
                            for (int i = imin + iloc * iinc;
                                 i < imin + (iloc + 1) * iinc; i += incx_)
                            {
                                buf3_ptr += dimz_;
                            }
                    }
                }
            }
            assert((buf3_ptr - &comm_buf3[0])
                   <= static_cast<int>(comm_buf3.size()));
        }
    }
    else
    { // grid_.mype_env().n_mpi_task(1)==1
        if (bc_[1] == 1)
        {

            // only for i already initialized
            for (int k = 0; k < nfunc_; k++)
            {
                ScalarType* pu = functions_[k]->uu();
                for (int j = 0; j < nghosts_ * incy_; j += incy_)
                    for (int i = nghosts_ * incx_;
                         i < (dimx_ + nghosts_) * incx_; i += incx_)
                    {
                        memcpy(&pu[i + nghosts_ + j],
                            &pu[i + nghosts_ + j + ymax], sdimz);
                        memcpy(&pu[i + nghosts_ * (incy_ + 1) + j + ymax],
                            &pu[i + nghosts_ * (incy_ + 1) + j], sdimz);
                    }
            }
        }
    }
    finishExchangeNorthSouth_tm_.stop();
}

template <typename ScalarType>
void GridFuncVector<ScalarType>::initiateUpDownComm(
    const int begin_color, const int end_color)
{
    std::vector<ScalarType>& comm_buf1(comm_buf1_[2]);
    std::vector<ScalarType>& comm_buf2(comm_buf2_[2]);
    std::vector<ScalarType>& comm_buf3(comm_buf3_[2]);
    std::vector<ScalarType>& comm_buf4(comm_buf4_[2]);

    const int ncolors = end_color - begin_color;

    const int sizebuffer = 1 + nfunc_max_global_ * up_down_size_;
    int iinit            = (dimy_ + 2 * nghosts_) * nghosts_;
    const int zmax       = dimz_;
    const int ione       = 1;

    if (down_)
    {
        grid_.mype_env().Irecv(
            &comm_buf3[0], sizebuffer, DOWN, &req_up_down_[1]);
    }
    if (up_)
    {
        grid_.mype_env().Irecv(&comm_buf4[0], sizebuffer, UP, &req_up_down_[3]);
    }

    if (up_)
    {
        ScalarType* buf1_ptr = &comm_buf1[0];
        // first element will tell how many functions (data) are in buffer
        *buf1_ptr = (ScalarType)ncolors;
        buf1_ptr++;
        for (int color = begin_color; color < end_color; color++)
        {
            for (short iloc = 0; iloc < nsubdivx_; iloc++)
            {
                const ScalarType* const uus
                    = functions_[color]->uu(nghosts_ + incy_ * iinit);
                *buf1_ptr = (ScalarType)gid_[iloc][color];
                buf1_ptr++;
                for (int j = 0; j < nghosts_; j++)
                {
                    Tcopy(&incxy_, &uus[zmax - 1 - j + iloc * incxy_ * incy_],
                        &incy_, buf1_ptr, &ione);
                    buf1_ptr += incxy_;
                }
            }
        }
        grid_.mype_env().Isend(
            &comm_buf1[0], 1 + ncolors * up_down_size_, UP, &req_up_down_[0]);
    }
    if (down_)
    {
        ScalarType* buf2_ptr = &comm_buf2[0];
        // first element will tell how many functions (data) are in buffer
        *buf2_ptr = (ScalarType)ncolors;
        buf2_ptr++;
        for (int color = begin_color; color < end_color; color++)
        {
            for (short iloc = 0; iloc < nsubdivx_; iloc++)
            {
                const ScalarType* const uus
                    = functions_[color]->uu(nghosts_ + incy_ * iinit);
                *buf2_ptr = (ScalarType)gid_[iloc][color];
                buf2_ptr++;
                for (int j = 0; j < nghosts_; j++)
                {
                    Tcopy(&incxy_, &uus[j + iloc * incxy_ * incy_], &incy_,
                        buf2_ptr, &ione);
                    buf2_ptr += incxy_;
                }
            }
        }
        grid_.mype_env().Isend(
            &comm_buf2[0], 1 + ncolors * up_down_size_, DOWN, &req_up_down_[2]);
    }
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::finishUpDownComm()
{
    finishExchangeUpDown_tm_.start();

    std::vector<ScalarType>& comm_buf3(comm_buf3_[2]);
    std::vector<ScalarType>& comm_buf4(comm_buf4_[2]);

    const int zmax = dimz_;

    if (grid_.mype_env().n_mpi_task(2) > 1)
    {

        int iinit      = (dimy_ + 2 * nghosts_) * nghosts_;
        const int ione = 1;
        if (down_)
        {
            ScalarType* buf3_ptr   = &comm_buf3[0];
            const int nremote_func = (int)(*buf3_ptr);
            buf3_ptr++;

            for (int k = 0; k < nremote_func; k++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const int gid = (int)(*buf3_ptr);
                    buf3_ptr++;
                    if (gid >= 0)
                    {
                        std::map<int, short>::const_iterator ilid
                            = gid2lid_.find(gid);
                        if (ilid != gid2lid_.end())
                        {
                            short lid = ilid->second;
                            if (gid == gid_[iloc][lid])
                            {
                                for (int j = 0; j < nghosts_; j++)
                                {
                                    ScalarType* const uus
                                        = functions_[lid]->uu(nghosts_ - 1 - j);
                                    Tcopy(&incxy_, buf3_ptr, &ione,
                                        &uus[incy_ * iinit
                                             + iloc * incxy_ * incy_],
                                        &incy_);
                                    buf3_ptr += incxy_;
                                }
                            }
                            else
                            {
                                buf3_ptr += incxy_ * nghosts_;
                            }
                        }
                        else
                        {
                            buf3_ptr += incxy_ * nghosts_;
                        }
                    }
                    else
                    {
                        buf3_ptr += incxy_ * nghosts_;
                    }
                }
            }
        }
        if (up_)
        {
            ScalarType* buf4_ptr   = &comm_buf4[0];
            const int nremote_func = (int)(*buf4_ptr);
            buf4_ptr++;

            for (int k = 0; k < nremote_func; k++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const int gid = (int)(*buf4_ptr);
                    buf4_ptr++;
                    if (gid >= 0)
                    {
                        std::map<int, short>::const_iterator ilid
                            = gid2lid_.find(gid);
                        if (ilid != gid2lid_.end())
                        {
                            short lid = ilid->second;
                            if (gid == gid_[iloc][lid])
                            {
                                for (int j = 0; j < nghosts_; j++)
                                {
                                    ScalarType* const uus = functions_[lid]->uu(
                                        nghosts_ + zmax + j);
                                    Tcopy(&incxy_, buf4_ptr, &ione,
                                        &uus[incy_ * iinit
                                             + iloc * incxy_ * incy_],
                                        &incy_);
                                    buf4_ptr += incxy_;
                                }
                            }
                            else
                            {
                                buf4_ptr += incxy_ * nghosts_;
                            }
                        }
                        else
                        {
                            buf4_ptr += incxy_ * nghosts_;
                        }
                    }
                    else
                    {
                        buf4_ptr += incxy_ * nghosts_;
                    }
                }
            }
        }
    }
    else
    {

        if (bc_[2] == 1) /* grid_.mype_env().n_mpi_task(2)==1 */
        {
            int iinit = (dimy_ + 2 * nghosts_) * nghosts_;

            for (int k = 0; k < nfunc_; k++)
            {
                ScalarType* pu = functions_[k]->uu();
                for (int j = 0; j < nghosts_; j++)
                {
                    Tcopy(&dimxy_, &pu[nghosts_ - j + iinit * incy_ + zmax - 1],
                        &incy_, &pu[nghosts_ - j + iinit * incy_ - 1], &incy_);
                }
                for (int j = 0; j < nghosts_; j++)
                {
                    Tcopy(&dimxy_, &pu[nghosts_ + j + iinit * incy_], &incy_,
                        &pu[nghosts_ + j + iinit * incy_ + zmax], &incy_);
                }
            }
        }
    }

    finishExchangeUpDown_tm_.stop();
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::initiateEastWestComm(
    const int begin_color, const int end_color)
{
    std::vector<ScalarType>& comm_buf1(comm_buf1_[0]);
    std::vector<ScalarType>& comm_buf2(comm_buf2_[0]);
    std::vector<ScalarType>& comm_buf3(comm_buf3_[0]);
    std::vector<ScalarType>& comm_buf4(comm_buf4_[0]);

    const int ncolors = end_color - begin_color;

    const int sizebuffer = 1 + nfunc_max_global_ * (east_west_size_ + 1);
    const int xmax       = dimx_ * grid_.inc(0);
    const size_t east_west_size_data = east_west_size_ * sizeof(ScalarType);

    /* Non-blocking MPI */
    if (east_)
        grid_.mype_env().Irecv(
            &comm_buf3[0], sizebuffer, EAST, &req_east_west_[3]);
    if (west_)
        grid_.mype_env().Irecv(
            &comm_buf4[0], sizebuffer, WEST, &req_east_west_[2]);
    if (west_)
    {
        ScalarType* buf1_ptr = &comm_buf1[0];
        // first element will tell how many functions (data) are in buffer
        *buf1_ptr = (ScalarType)ncolors;
        buf1_ptr++;

        const int initu = east_west_size_;
        for (int color = begin_color; color < end_color; color++)
        {
            *buf1_ptr = (ScalarType)gid_[0][color];
            buf1_ptr++;
            const ScalarType* const pu = functions_[color]->uu(initu);
            memcpy(buf1_ptr, pu, east_west_size_data);
            buf1_ptr += east_west_size_;
        }
        grid_.mype_env().Isend(
            &comm_buf1[0], sizebuffer, WEST, &req_east_west_[1]);
    }
    if (east_)
    {
        ScalarType* buf2_ptr = &comm_buf2[0];
        // first element will tell how many functions (data) are in buffer
        *buf2_ptr = (ScalarType)ncolors;
        buf2_ptr++;

        const int initu = xmax;
        for (int color = begin_color; color < end_color; color++)
        {
            *buf2_ptr = (ScalarType)gid_[nsubdivx_ - 1][color];
            buf2_ptr++;
            const ScalarType* const pu = functions_[color]->uu(initu);
            memcpy(buf2_ptr, pu, east_west_size_data);
            buf2_ptr += east_west_size_;
        }
        grid_.mype_env().Isend(
            &comm_buf2[0], sizebuffer, EAST, &req_east_west_[0]);
    }
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::finishEastWestComm()
{
    finishExchangeEastWest_tm_.start();

    std::vector<ScalarType>& comm_buf3(comm_buf3_[0]);
    std::vector<ScalarType>& comm_buf4(comm_buf4_[0]);

    const int xmax = dimx_ * grid_.inc(0);

    if (grid_.mype_env().n_mpi_task(0) > 1)
    {

        const size_t east_west_size_data = east_west_size_ * sizeof(ScalarType);
        if (east_)
        {
            ScalarType* buf3_ptr   = &comm_buf3[0];
            const int nremote_func = (int)(*buf3_ptr);
            buf3_ptr++;

            const int initu = xmax + east_west_size_;
            for (int k = 0; k < nremote_func; k++)
            {
                int gid = (int)(*buf3_ptr);
                buf3_ptr++;
                if (gid >= 0)
                {
                    std::map<int, short>::const_iterator ilid
                        = gid2lid_.find(gid);
                    if (ilid != gid2lid_.end())
                    {
                        short lid = ilid->second;
                        if (gid == gid_[nsubdivx_ - 1][lid])
                        {
                            ScalarType* const pu = functions_[lid]->uu(initu);
                            memcpy(pu, buf3_ptr, east_west_size_data);
                        }
                    }
                }
                buf3_ptr += east_west_size_;
            }
        }
        if (west_)
        {
            ScalarType* buf4_ptr   = &comm_buf4[0];
            const int nremote_func = (int)(*buf4_ptr);
            buf4_ptr++;

            const int initu = 0;
            for (int k = 0; k < nremote_func; k++)
            {
                int gid = (int)(*buf4_ptr);
                buf4_ptr++;
                if (gid >= 0)
                {
                    std::map<int, short>::const_iterator ilid
                        = gid2lid_.find(gid);
                    if (ilid != gid2lid_.end())
                    {
                        short lid = ilid->second;
                        if (gid == gid_[0][lid])
                        {
                            ScalarType* const pu = functions_[lid]->uu(initu);
                            memcpy(pu, buf4_ptr, east_west_size_data);
                        }
                    }
                }
                buf4_ptr += east_west_size_;
            }
        }
    }
    else
    { /* grid_.mype_env().n_mpi_task(0)==1 */

        if (bc_[0] == 1)
        {
            const size_t east_west_size_data
                = east_west_size_ * sizeof(ScalarType);

            for (int k = 0; k < nfunc_; k++)
            {
                ScalarType* pu = functions_[k]->uu();
                memcpy(&pu[0], &pu[xmax], east_west_size_data);
                memcpy(&pu[east_west_size_ + xmax], &pu[east_west_size_],
                    east_west_size_data);
            }
        }
    }
    finishExchangeEastWest_tm_.stop();
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::trade_boundaries()
{
    if (updated_boundaries_) return;

    for (int k = 0; k < nfunc_; k++)
    {
        assert(functions_[k]->updated_boundaries()
               == functions_[0]->updated_boundaries());
        assert(functions_[k]->updated_boundaries() == updated_boundaries_);
    }
    assert(dimx_ >= nghosts_);
    assert(dimy_ >= nghosts_);
    assert(dimz_ >= nghosts_);

    trade_bc_tm_.start();

    bool direction[3] = { !bc_[0], !bc_[1], !bc_[2] };
    if (bc_[0] != 1 || bc_[1] != 1 || bc_[2] != 1)
        for (int i = 0; i < nfunc_; i++)
            functions_[i]->setBoundaryValues(0., direction);

#if 0
    for(int k=0;k<nfunc_;k+=5)
       trade_boundaries_colors(k,k+5);
    
    grid_.mype_env().barrier(); 
    
    return;
#endif

    // Create buffers
    allocate_buffers(nfunc_max_global_);

    if (grid_.mype_env().n_mpi_task(1) > 1)
    {
        initiateNorthSouthComm(0, nfunc_);
    }
    if (!skinny_stencil_)
    {
        wait_north_south();
        finishNorthSouthComm();
    }

    if (grid_.mype_env().n_mpi_task(2) > 1)
    {
        initiateUpDownComm(0, nfunc_);
    }
    if (!skinny_stencil_)
    {
        wait_up_down();
        finishUpDownComm();
    }

    if (grid_.mype_env().n_mpi_task(0) > 1)
    {
        initiateEastWestComm(0, nfunc_);
    }

    if (!skinny_stencil_)
    {
        wait_east_west();
        finishEastWestComm();
    }

    if (skinny_stencil_)
    {
        wait_north_south();
        finishNorthSouthComm();

        wait_up_down();
        finishUpDownComm();

        wait_east_west();
        finishEastWestComm();
    }

    updated_boundaries_ = true;

    for (int k = 0; k < nfunc_; k++)
        functions_[k]->set_updated_boundaries(true);

    trade_bc_tm_.stop();
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::restrict3D(GridFuncVector& ucoarse)
{
    if (!updated_boundaries_) trade_boundaries();

    for (short k = 0; k < nfunc_; k++)
    {
        functions_[k]->restrict3D(*ucoarse.functions_[k]);
    }
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::extend3D(GridFuncVector& ucoarse)
{
    if (!ucoarse.updated_boundaries_) ucoarse.trade_boundaries();

    for (short k = 0; k < nfunc_; k++)
    {
        functions_[k]->extend3D(*ucoarse.functions_[k]);
    }

    updated_boundaries_ = false;
}
template <typename ScalarType>
GridFuncVector<ScalarType>& GridFuncVector<ScalarType>::operator-=(
    const GridFuncVector& func)
{
    assert(func.grid_.sizeg() == grid_.sizeg());
    assert(func.grid_.ghost_pt() == grid_.ghost_pt());
    assert(this != &func);

    for (short k = 0; k < nfunc_; k++)
    {
        (*functions_[k]) -= (*func.functions_[k]);
    }

    updated_boundaries_ = (func.updated_boundaries_ && updated_boundaries_);

    return *this;
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::axpy(
    const double alpha, const GridFuncVector<ScalarType>& func)
{
    for (short k = 0; k < nfunc_; k++)
    {
        functions_[k]->axpy(alpha, *func.functions_[k]);
    }

    updated_boundaries_ = (func.updated_boundaries_ && updated_boundaries_);
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::init_vect(
    const int k, ScalarType* vv, const char dis) const
{
    functions_[k]->init_vect(vv, dis);
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::getValues(const int k, float* vv) const
{
    assert(k < static_cast<int>(functions_.size()));
    functions_[k]->getValues(vv);
}
template <typename ScalarType>
void GridFuncVector<ScalarType>::getValues(const int k, double* vv) const
{
    assert(k < static_cast<int>(functions_.size()));
    functions_[k]->getValues(vv);
}

// build list of local gids I need ghost values for
template <typename ScalarType>
void GridFuncVector<ScalarType>::communicateRemoteGids(
    const int begin_color, const int end_color)
{
    const int ncolors = end_color - begin_color;

    MPI_Request gid_req[6];

    // build list of local gids I need ghost values for
    std::vector<int> local_gids;
    const short ndata = ncolors * nsubdivx_;
    local_gids.reserve(ndata);
    for (short color = begin_color; color < end_color; color++)
    {
        for (short iloc = 0; iloc < nsubdivx_; iloc++)
        {
            if (color < static_cast<int>(gid_[iloc].size()))
                local_gids.push_back(gid_[iloc][color]);
            else
                local_gids.push_back(-1);
        }
    }

    for (short i = 0; i < 6; i++)
        remote_gids_[i].resize(ndata);

    if (north_)
        grid_.mype_env().Irecv(
            &remote_gids_[NORTH][0], ndata, NORTH, &gid_req[NORTH]);
    if (south_)
        grid_.mype_env().Irecv(
            &remote_gids_[SOUTH][0], ndata, SOUTH, &gid_req[SOUTH]);
    if (east_)
        grid_.mype_env().Irecv(
            &remote_gids_[EAST][0], ndata, EAST, &gid_req[EAST]);
    if (west_)
        grid_.mype_env().Irecv(
            &remote_gids_[WEST][0], ndata, WEST, &gid_req[WEST]);
    if (up_)
        grid_.mype_env().Irecv(&remote_gids_[UP][0], ndata, UP, &gid_req[UP]);
    if (down_)
        grid_.mype_env().Irecv(
            &remote_gids_[DOWN][0], ndata, DOWN, &gid_req[DOWN]);

    if (north_)
        grid_.mype_env().Isend(&local_gids[0], ndata, NORTH, &gid_req[NORTH]);
    if (south_)
        grid_.mype_env().Isend(&local_gids[0], ndata, SOUTH, &gid_req[SOUTH]);
    if (east_)
        grid_.mype_env().Isend(&local_gids[0], ndata, EAST, &gid_req[EAST]);
    if (west_)
        grid_.mype_env().Isend(&local_gids[0], ndata, WEST, &gid_req[WEST]);
    if (up_) grid_.mype_env().Isend(&local_gids[0], ndata, UP, &gid_req[UP]);
    if (down_)
        grid_.mype_env().Isend(&local_gids[0], ndata, DOWN, &gid_req[DOWN]);

    if (south_) MPI_Wait(gid_req + SOUTH, MPI_STATUS_IGNORE);
    if (north_) MPI_Wait(gid_req + NORTH, MPI_STATUS_IGNORE);
    if (west_) MPI_Wait(gid_req + WEST, MPI_STATUS_IGNORE);
    if (east_) MPI_Wait(gid_req + EAST, MPI_STATUS_IGNORE);
    if (down_) MPI_Wait(gid_req + DOWN, MPI_STATUS_IGNORE);
    if (up_) MPI_Wait(gid_req + UP, MPI_STATUS_IGNORE);
}

// assumes all processors call with same arguments
template <typename ScalarType>
void GridFuncVector<ScalarType>::trade_boundaries_colors(
    const short first_color, const short last_color)
{
    assert(first_color >= 0);
    assert(last_color >= 0);
    assert(comm_buf3_.size() > 0);

    grid_.mype_env().barrier();
    // if (onpe0)
    //    std::cout << "Color " << first_color << " to " << last_color - 1
    //              << std::endl;

    if (updated_boundaries_) return;

    std::vector<ScalarType>& comm_buf1(comm_buf1_[0]);
    std::vector<ScalarType>& comm_buf2(comm_buf2_[0]);
    std::vector<ScalarType>& comm_buf3(comm_buf3_[0]);
    std::vector<ScalarType>& comm_buf4(comm_buf4_[0]);

    const short ncolors   = last_color - first_color;
    const short end_color = std::min(last_color, (short)functions_.size());

    for (short color = first_color + 1; color < end_color; color++)
        assert(functions_[color]->updated_boundaries()
               == functions_[first_color]->updated_boundaries());

    trade_bc_colors_tm_.start();

    const int xmax = dimx_ * grid_.inc(0);
    const int ymax = dimy_ * grid_.inc(1);
    const int zmax = dimz_;

    const size_t sdimz = dimz_ * sizeof(ScalarType);

    communicateRemoteGids(first_color, last_color);

    const int ione = 1;

    // Create buffers
    allocate_buffers(ncolors);

    if (grid_.mype_env().n_mpi_task(1) > 1)
    {

        const size_t sizebuffer = 1 + ncolors * south_north_size_;
        if (north_)
            grid_.mype_env().Irecv(
                &comm_buf4[0], sizebuffer, NORTH, &req_north_south_[3]);
        if (south_)
            grid_.mype_env().Irecv(
                &comm_buf3[0], sizebuffer, SOUTH, &req_north_south_[1]);

        const int imin = nghosts_ * incx_;
        const int iinc = incx_ * dimx_ / nsubdivx_;

        const int jmax = nghosts_ * incy_;

        if (south_)
        {
            ScalarType* buf2_ptr = &comm_buf2[0];
            // first element will tell how many functions (data) are in buffer
            *buf2_ptr = (ScalarType)ncolors;
            buf2_ptr++;
            // pack data
            for (short color = first_color; color < last_color; color++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const short rgid
                        = remote_gids_[SOUTH][nsubdivx_ * color + iloc];
                    const short k = gid2lid_[rgid];
                    assert(k < static_cast<int>(functions_.size()));
                    const ScalarType* uus
                        = functions_[k]->uu(nghosts_ * (incy_ + 1));
                    *buf2_ptr = (ScalarType)gid_[iloc][k];
                    buf2_ptr++;
                    for (int j = 0; j < jmax; j += incy_)
                        for (int i = imin + iloc * iinc;
                             i < imin + (iloc + 1) * iinc; i += incx_)
                        {
                            memcpy(buf2_ptr, &uus[i + j], sdimz);
                            buf2_ptr += dimz_;
                        }
                }
            }
            grid_.mype_env().Isend(&comm_buf2[0],
                1 + ncolors * south_north_size_, SOUTH, &req_north_south_[2]);
        }

        if (north_)
        {
            ScalarType* buf1_ptr = &comm_buf1[0];
            // first element will tell how many functions (data) are in buffer
            *buf1_ptr = (ScalarType)ncolors;
            buf1_ptr++;
            // pack data
            for (short color = first_color; color < last_color; color++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const short rgid
                        = remote_gids_[NORTH][nsubdivx_ * color + iloc];
                    const short k         = gid2lid_[rgid];
                    const ScalarType* uus = functions_[k]->uu(nghosts_ + ymax);
                    *buf1_ptr             = (ScalarType)gid_[iloc][k];
                    buf1_ptr++;
                    for (int j = 0; j < jmax; j += incy_)
                        for (int i = imin + iloc * iinc;
                             i < imin + (iloc + 1) * iinc; i += incx_)
                        {
                            memcpy(buf1_ptr, &uus[i + j], sdimz);
                            buf1_ptr += dimz_;
                        }
                }
            }
            grid_.mype_env().Isend(&comm_buf1[0],
                1 + ncolors * south_north_size_, NORTH, &req_north_south_[0]);
        }

        wait_north_south();

        if (north_)
        {
            ScalarType* buf4_ptr   = &comm_buf4[0];
            const int nremote_func = (int)(*buf4_ptr);
            buf4_ptr++;

            for (int k = 0; k < nremote_func; k++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const int gid = (int)(*buf4_ptr);
                    buf4_ptr++;

                    if (gid >= 0)
                    {
                        std::map<int, short>::const_iterator ilid
                            = gid2lid_.find(gid);
                        if (ilid != gid2lid_.end())
                        {
                            const short lid = ilid->second;
                            if (gid == gid_[iloc][lid])
                            {
                                ScalarType* uus = functions_[lid]->uu(
                                    nghosts_ * (incy_ + 1) + ymax);
                                for (int j = 0; j < jmax; j += incy_)
                                    for (int i = imin + iloc * iinc;
                                         i < imin + (iloc + 1) * iinc;
                                         i += incx_)
                                    {
                                        memcpy(&uus[i + j], buf4_ptr, sdimz);
                                        buf4_ptr += dimz_;
                                    }
                            }
                            else
                            {
                                for (int j = 0; j < jmax; j += incy_)
                                    for (int i = imin + iloc * iinc;
                                         i < imin + (iloc + 1) * iinc;
                                         i += incx_)
                                        buf4_ptr += dimz_;
                            }
                        }
                        else
                        {
                            for (int j = 0; j < jmax; j += incy_)
                                for (int i = imin + iloc * iinc;
                                     i < imin + (iloc + 1) * iinc; i += incx_)
                                    buf4_ptr += dimz_;
                        }
                    }
                    else
                    {
                        for (int j = 0; j < jmax; j += incy_)
                            for (int i = imin + iloc * iinc;
                                 i < imin + (iloc + 1) * iinc; i += incx_)
                            {
                                buf4_ptr += dimz_;
                            }
                    }
                }
            }
            assert((buf4_ptr - &comm_buf4[0]) <= static_cast<int>(sizebuffer));
        }

        if (south_)
        {
            ScalarType* buf3_ptr   = &comm_buf3[0];
            const int nremote_func = (int)(*buf3_ptr);
            buf3_ptr++;

            for (int k = 0; k < nremote_func; k++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const int gid = (int)(*buf3_ptr);
                    buf3_ptr++;
                    if (gid >= 0)
                    {
                        std::map<int, short>::const_iterator ilid
                            = gid2lid_.find(gid);
                        if (ilid != gid2lid_.end())
                        {
                            const short lid = ilid->second;
                            if (gid == gid_[iloc][lid])
                            {
                                ScalarType* uus = functions_[lid]->uu(nghosts_);
                                for (int j = 0; j < jmax; j += incy_)
                                    for (int i = imin + iloc * iinc;
                                         i < imin + (iloc + 1) * iinc;
                                         i += incx_)
                                    {
                                        memcpy(&uus[i + j], buf3_ptr, sdimz);
                                        buf3_ptr += dimz_;
                                    }
                            }
                            else
                            {
                                for (int j = 0; j < jmax; j += incy_)
                                    for (int i = imin + iloc * iinc;
                                         i < imin + (iloc + 1) * iinc;
                                         i += incx_)
                                        buf3_ptr += dimz_;
                            }
                        }
                        else
                        {
                            for (int j = 0; j < jmax; j += incy_)
                                for (int i = imin + iloc * iinc;
                                     i < imin + (iloc + 1) * iinc; i += incx_)
                                    buf3_ptr += dimz_;
                        }
                    }
                    else
                    {
                        for (int j = 0; j < jmax; j += incy_)
                            for (int i = imin + iloc * iinc;
                                 i < imin + (iloc + 1) * iinc; i += incx_)
                            {
                                buf3_ptr += dimz_;
                            }
                    }
                }
            }
            assert((buf3_ptr - &comm_buf3[0]) <= static_cast<int>(sizebuffer));
        }
    }
    else if (bc_[1] == 1)
    { // grid_.mype_env().n_mpi_task(1)==1

        // only for i already initialized
        for (short color = first_color; color < end_color; color++)
        {
            ScalarType* pu = functions_[color]->uu();
            for (int j = 0; j < nghosts_ * incy_; j += incy_)
                for (int i = nghosts_ * incx_; i < (dimx_ + nghosts_) * incx_;
                     i += incx_)
                {

                    memcpy(&pu[i + nghosts_ + j], &pu[i + nghosts_ + j + ymax],
                        sdimz);
                    memcpy(&pu[i + nghosts_ * (incy_ + 1) + j + ymax],
                        &pu[i + nghosts_ * (incy_ + 1) + j], sdimz);
                }
        }
    }

    int iinit = (dimy_ + 2 * nghosts_) * nghosts_;

    if (grid_.mype_env().n_mpi_task(2) > 1)
    {

        const int sizebuffer = 1 + ncolors * up_down_size_;
        if (down_)
        {
            assert(sizebuffer <= static_cast<int>(comm_buf3.size()));
            grid_.mype_env().Irecv(
                &comm_buf3[0], sizebuffer, DOWN, &req_up_down_[1]);
        }
        if (up_)
        {
            grid_.mype_env().Irecv(
                &comm_buf4[0], sizebuffer, UP, &req_up_down_[3]);
        }

        if (up_)
        {
            ScalarType* buf1_ptr = &comm_buf1[0];
            // first element will tell how many functions (data) are in buffer
            *buf1_ptr = (ScalarType)ncolors;
            buf1_ptr++;
            for (short color = first_color; color < last_color; color++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const short rgid
                        = remote_gids_[UP][nsubdivx_ * color + iloc];
                    const short k = gid2lid_[rgid];
                    const ScalarType* const uus
                        = functions_[k]->uu(nghosts_ + incy_ * iinit);
                    *buf1_ptr = (ScalarType)gid_[iloc][k];
                    buf1_ptr++;
                    for (int j = 0; j < nghosts_; j++)
                    {
                        Tcopy(&incxy_,
                            &uus[zmax - 1 - j + iloc * incxy_ * incy_], &incy_,
                            buf1_ptr, &ione);
                        buf1_ptr += incxy_;
                    }
                }
            }
            grid_.mype_env().Isend(&comm_buf1[0], 1 + ncolors * up_down_size_,
                UP, &req_up_down_[0]);
        }
        if (down_)
        {
            ScalarType* buf2_ptr = &comm_buf2[0];
            // first element will tell how many functions (data) are in buffer
            *buf2_ptr = (ScalarType)ncolors;
            buf2_ptr++;
            for (short color = first_color; color < last_color; color++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const short rgid
                        = remote_gids_[DOWN][nsubdivx_ * color + iloc];
                    const short k = gid2lid_[rgid];
                    const ScalarType* const uus
                        = functions_[k]->uu(nghosts_ + incy_ * iinit);
                    *buf2_ptr = (ScalarType)gid_[iloc][k];
                    buf2_ptr++;
                    for (int j = 0; j < nghosts_; j++)
                    {
                        Tcopy(&incxy_, &uus[j + iloc * incxy_ * incy_], &incy_,
                            buf2_ptr, &ione);
                        buf2_ptr += incxy_;
                    }
                }
            }
            grid_.mype_env().Isend(&comm_buf2[0], 1 + ncolors * up_down_size_,
                DOWN, &req_up_down_[2]);
        }

        wait_up_down();

        if (down_)
        {
            ScalarType* buf3_ptr   = &comm_buf3[0];
            const int nremote_func = (int)(*buf3_ptr);
            buf3_ptr++;

            for (int k = 0; k < nremote_func; k++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const int gid = (int)(*buf3_ptr);
                    buf3_ptr++;
                    if (gid >= 0)
                    {
                        std::map<int, short>::const_iterator ilid
                            = gid2lid_.find(gid);
                        if (ilid != gid2lid_.end())
                        {
                            short lid = ilid->second;
                            if (gid == gid_[iloc][lid])
                            {
                                for (int j = 0; j < nghosts_; j++)
                                {
                                    ScalarType* const uus
                                        = functions_[lid]->uu(nghosts_ - 1 - j);
                                    Tcopy(&incxy_, buf3_ptr, &ione,
                                        &uus[incy_ * iinit
                                             + iloc * incxy_ * incy_],
                                        &incy_);
                                    buf3_ptr += incxy_;
                                }
                            }
                            else
                            {
                                buf3_ptr += incxy_ * nghosts_;
                            }
                        }
                        else
                        {
                            buf3_ptr += incxy_ * nghosts_;
                        }
                    }
                    else
                    {
                        buf3_ptr += incxy_ * nghosts_;
                    }
                }
            }
        }
        if (up_)
        {
            ScalarType* buf4_ptr   = &comm_buf4[0];
            const int nremote_func = (int)(*buf4_ptr);
            buf4_ptr++;

            for (int k = 0; k < nremote_func; k++)
            {
                for (short iloc = 0; iloc < nsubdivx_; iloc++)
                {
                    const int gid = (int)(*buf4_ptr);
                    buf4_ptr++;
                    if (gid >= 0)
                    {
                        std::map<int, short>::const_iterator ilid
                            = gid2lid_.find(gid);
                        if (ilid != gid2lid_.end())
                        {
                            short lid = ilid->second;
                            if (gid == gid_[iloc][lid])
                            {
                                for (int j = 0; j < nghosts_; j++)
                                {
                                    ScalarType* const uus = functions_[lid]->uu(
                                        nghosts_ + zmax + j);
                                    Tcopy(&incxy_, buf4_ptr, &ione,
                                        &uus[incy_ * iinit
                                             + iloc * incxy_ * incy_],
                                        &incy_);
                                    buf4_ptr += incxy_;
                                }
                            }
                            else
                            {
                                buf4_ptr += incxy_ * nghosts_;
                            }
                        }
                        else
                        {
                            buf4_ptr += incxy_ * nghosts_;
                        }
                    }
                    else
                    {
                        buf4_ptr += incxy_ * nghosts_;
                    }
                }
            }
        }
    }
    else if (bc_[2] == 1)
    { /* grid_.mype_env().n_mpi_task(2)==1 */

        for (short color = first_color; color < end_color; color++)
        {
            ScalarType* pu = functions_[color]->uu();
            for (int j = 0; j < nghosts_; j++)
            {
                Tcopy(&dimxy_, &pu[nghosts_ - j + iinit * incy_ + zmax - 1],
                    &incy_, &pu[nghosts_ - j + iinit * incy_ - 1], &incy_);
            }
            for (int j = 0; j < nghosts_; j++)
            {
                Tcopy(&dimxy_, &pu[nghosts_ + j + iinit * incy_], &incy_,
                    &pu[nghosts_ + j + iinit * incy_ + zmax], &incy_);
            }
        }
    }

    const int east_west_size_        = nghosts_ * incx_;
    const size_t east_west_size_data = east_west_size_ * sizeof(ScalarType);

    if (grid_.mype_env().n_mpi_task(0) > 1)
    {
        const int sizebuffer = 1 + ncolors * (east_west_size_ + 1);

        /* Non-blocking MPI */
        if (east_)
            grid_.mype_env().Irecv(
                &comm_buf3[0], sizebuffer, EAST, &req_east_west_[3]);
        if (west_)
            grid_.mype_env().Irecv(
                &comm_buf4[0], sizebuffer, WEST, &req_east_west_[2]);

        if (west_)
        {
            ScalarType* buf1_ptr = &comm_buf1[0];
            // first element will tell how many functions (data) are in buffer
            *buf1_ptr = (ScalarType)ncolors;
            buf1_ptr++;

            const int initu = east_west_size_;
            for (short color = first_color; color < last_color; color++)
            {
                const short rgid
                    = remote_gids_[WEST][nsubdivx_ * color + nsubdivx_ - 1];
                if (rgid >= 0)
                {
                    const short k = gid2lid_[rgid];
                    *buf1_ptr     = (ScalarType)gid_[0][k];
                    buf1_ptr++;
                    if (k < static_cast<int>(functions_.size()))
                    {
                        const ScalarType* const pu = functions_[k]->uu(initu);
                        memcpy(buf1_ptr, pu, east_west_size_data);
                    }
                }
                else
                {
                    *buf1_ptr = -1.;
                    buf1_ptr++;
                }
                buf1_ptr += east_west_size_;
            }
            grid_.mype_env().Isend(
                &comm_buf1[0], sizebuffer, WEST, &req_east_west_[1]);
        }
        if (east_)
        {
            ScalarType* buf2_ptr = &comm_buf2[0];
            // first element will tell how many functions (data) are in buffer
            *buf2_ptr = (ScalarType)ncolors;
            buf2_ptr++;

            const int initu = xmax;
            for (short color = first_color; color < last_color; color++)
            {
                const short rgid = remote_gids_[EAST][nsubdivx_ * color + 0];
                if (rgid >= 0)
                {
                    const short k = gid2lid_[rgid];
                    *buf2_ptr     = (ScalarType)gid_[nsubdivx_ - 1][k];
                    buf2_ptr++;
                    if (k < static_cast<int>(functions_.size()))
                    {
                        const ScalarType* const pu = functions_[k]->uu(initu);
                        memcpy(buf2_ptr, pu, east_west_size_data);
                    }
                }
                else
                {
                    *buf2_ptr = -1.;
                    buf2_ptr++;
                }
                buf2_ptr += east_west_size_;
            }
            grid_.mype_env().Isend(
                &comm_buf2[0], sizebuffer, EAST, &req_east_west_[0]);
        }

        wait_east_west();

        if (east_)
        {
            ScalarType* buf3_ptr   = &comm_buf3[0];
            const int nremote_func = (int)(*buf3_ptr);
            buf3_ptr++;

            const int initu = xmax + east_west_size_;
            for (int k = 0; k < nremote_func; k++)
            {
                int gid = (int)(*buf3_ptr);
                buf3_ptr++;
                if (gid >= 0)
                {
                    std::map<int, short>::const_iterator ilid
                        = gid2lid_.find(gid);
                    if (ilid != gid2lid_.end())
                    {
                        short lid = ilid->second;
                        if (gid == gid_[nsubdivx_ - 1][lid])
                        {
                            ScalarType* const pu = functions_[lid]->uu(initu);
                            memcpy(pu, buf3_ptr, east_west_size_data);
                        }
                    }
                }
                buf3_ptr += east_west_size_;
            }
        }
        if (west_)
        {
            ScalarType* buf4_ptr   = &comm_buf4[0];
            const int nremote_func = (int)(*buf4_ptr);
            buf4_ptr++;

            const int initu = 0;
            for (int k = 0; k < nremote_func; k++)
            {
                int gid = (int)(*buf4_ptr);
                buf4_ptr++;
                if (gid >= 0)
                {
                    std::map<int, short>::const_iterator ilid
                        = gid2lid_.find(gid);
                    if (ilid != gid2lid_.end())
                    {
                        short lid = ilid->second;
                        if (gid == gid_[0][lid])
                        {
                            ScalarType* const pu = functions_[lid]->uu(initu);
                            memcpy(pu, buf4_ptr, east_west_size_data);
                        }
                    }
                }
                buf4_ptr += east_west_size_;
            }
        }
    }
    else if (bc_[0] == 1)
    { /* grid_.mype_env().n_mpi_task(0)==1 */

        for (short color = first_color; color < end_color; color++)
        {
            ScalarType* pu = functions_[color]->uu();
            memcpy(&pu[0], &pu[xmax], east_west_size_data);
            memcpy(&pu[east_west_size_ + xmax], &pu[east_west_size_],
                east_west_size_data);
        }
    }

    updated_boundaries_ = true;

    for (short color = first_color; color < end_color; color++)
        functions_[color]->set_updated_boundaries(true);

    trade_bc_colors_tm_.stop();
}

template class GridFuncVector<double>;
template class GridFuncVector<float>;
template void GridFuncVector<float>::assign(
    const int i, const float* const v, const char dis);
template void GridFuncVector<float>::assign(
    const int i, const double* const v, const char dis);
template void GridFuncVector<double>::assign(
    const int i, const float* const v, const char dis);
template void GridFuncVector<double>::assign(
    const int i, const double* const v, const char dis);

} // namespace pb
