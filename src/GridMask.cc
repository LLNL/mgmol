// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cassert>
#include <iostream>

#include "Control.h"
#include "GridFunc.h"
#include "GridMask.h"
#include "SubCell.h"

// #define DEBUG 1

Timer GridMask::init_tm_("GridMask::init");

std::vector<SubCell*> GridMask::sub_cell_;
int GridMask::ninstances_ = 0;

// Constructor
GridMask::GridMask(const unsigned short nclevels, const unsigned short subdivx,
    const pb::Grid& mygrid)
    : grid_(mygrid), subdivx_(subdivx)
{

    int numpt = grid_.size();
    assert(numpt > 0);
    assert(numpt < 1000000000);
    assert(subdivx > 0);
    assert(subdivx < 100);
    assert(nclevels < 10);

    nlevels_ = nclevels + 1;
    if (nlevels_ == 0) nlevels_ = 1;
    const double h[3] = { grid_.hgrid(0), grid_.hgrid(1), grid_.hgrid(2) };
    // width of buffer region between localized orbitals
    delta_ = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]) + 1.e-8;
    // Take into account Differential operator
    delta_ *= (double)grid_.ghost_pt();

    assert(nlevels_ > 0);

    loc_numpt_.resize(nlevels_);
    for (unsigned short l = 0; l < nlevels_; l++)
    {
        loc_numpt_[l] = (numpt >> (3 * l)) / subdivx;
        assert(loc_numpt_[l]
               == numpt / (subdivx * (1 << l) * (1 << l) * (1 << l)));
    }
    subdim0_.resize(nlevels_);
    for (unsigned short l = 0; l < nlevels_; l++)
    {
        subdim0_[l] = (grid_.dim(0) >> l) / subdivx_;
    }
    mask_not_zero_.resize(nlevels_);
    for (unsigned short l = 0; l < nlevels_; l++)
    {
        mask_not_zero_[l].resize(subdivx_);
        // default value: all the masks are at 1
        for (unsigned short i = 0; i < subdivx_; i++)
        {
            mask_not_zero_[l][i] = 1;
        }
    }
    lmask_.resize(nlevels_);
    for (unsigned short l = 0; l < nlevels_; l++)
    {
        lmask_[l].resize(subdivx_);
    }

    // allocate static data
    if (sub_cell_.empty())
    {
        sub_cell_.resize(nlevels_);
        for (unsigned short l = 0; l < nlevels_; l++)
            sub_cell_[l] = new SubCell(grid_, subdivx, l);
    }

    ninstances_++;
}

// Destructor
GridMask::~GridMask()
{
    if (ninstances_ == 1)
    {
        for (unsigned short i = 0; i < sub_cell_.size(); i++)
            delete sub_cell_[i];
    }
    ninstances_--;
}

GridMask& GridMask::operator=(const GridMask& grid_mask)
{
    if (this == &grid_mask) return *this;

    nlevels_   = grid_mask.nlevels_;
    subdivx_   = grid_mask.subdivx_;
    delta_     = grid_mask.delta_;
    loc_numpt_ = grid_mask.loc_numpt_;
    subdim0_   = grid_mask.subdim0_;

    center_ = grid_mask.center_;
    radius_ = grid_mask.radius_;

    lmask_ = grid_mask.lmask_;

    mask_not_zero_ = grid_mask.mask_not_zero_;

    return *this;
}

// Initialize mask of radius "rcut" of center "v",
// at level "level"
int GridMask::init(const Vector3D& center, const double rcut,
    const unsigned short level, lmasktype (*func)(const double))
{
    assert(subdivx_ > 0);

    int n = 0;
    for (unsigned short iloc = 0; iloc < subdivx_; iloc++)
        n += init(center, rcut, iloc, level, func);
    return n;
}

// Initialize mask of radius "rcut" of center "center" in region "iloc",
// at level "level"
int GridMask::init(const Vector3D& center, const double rcut,
    const unsigned short iloc, const unsigned short level,
    lmasktype (*func)(const double))
{
    assert(rcut > 0.001);
    assert((*func)(0) == 1.);
    assert(subdim0_[level] > 0);

#pragma omp master
    init_tm_.start();

    //(*MPIdata::sout)<<"GridMask::init(), center="<<center<<",
    // rcut="<<rcut<<endl;
    Control& ct = *(Control::instance());
    center_     = center;
    radius_     = rcut;

    assert(subdim0_[level] > 1);
    const int dim1 = (grid_.dim(1) >> level);
    const int dim2 = (grid_.dim(2) >> level);
    const int size = loc_numpt_[level];

    std::vector<lmasktype> mask(size);

    double delta = delta_;
    for (unsigned short l = 0; l < level; l++)
        delta *= 2.;

    double overlap_radius = rcut + delta;
    bool overlap          = sub_cell_[level]->spherePossibleOverlap(
        center, overlap_radius, iloc, ct.bcPoisson);

    // try cheap cases first
    if (!overlap)
    { // Mask=0 in whole subdomain
        assign(iloc, level, 0, 0, mask);
        init_tm_.stop();
        return 0;
    }
    else if (100. < rcut)
    { // Mask =1
        assign1(iloc, level);
        init_tm_.stop();
        return size;
    }

    const double rc_inv = 1. / rcut;
    const Vector3D ll(grid_.ll(0), grid_.ll(1), grid_.ll(2));

    const double start[3] = { grid_.start(0), grid_.start(1), grid_.start(2) };

    const double hgrid[3] = { grid_.hgrid(0) * (1 << level),
        grid_.hgrid(1) * (1 << level), grid_.hgrid(2) * (1 << level) };

#ifdef DEBUG
    int incy   = dim2;
    int incx   = dim2 * dim1;
    int lcount = 0;
    //(*MPIdata::sout)<<"dim=("<<subdim0_[level]<<", "<<dim1<<",
    //"<<dim2<<")"<<endl;
#endif
    int icount = 0; // nb. grid points inside radius overlap_radius
    int ncount = 0; // nb. grid points inside radius rcut
    Vector3D xc;
    xc[0]            = start[0] + iloc * subdim0_[level] * hgrid[0];
    lmasktype* pmask = &mask[0];
    for (int ix = 0; ix < subdim0_[level]; ix++)
    {
        xc[1] = start[1];

        for (int iy = 0; iy < dim1; iy++)
        {
            xc[2] = start[2];

            for (int iz = 0; iz < dim2; iz++)
            {
                assert(sub_cell_[level]->includes(xc, iloc));

                const double r = xc.minimage(center, ll, ct.bcPoisson);
                assert(r >= 0.);
                assert(r * r
                       < ll[0] * ll[0] + ll[1] * ll[1] + ll[2] * ll[2] + 0.1);
                if (r < overlap_radius)
                {
                    icount++;

                    if (r <= rcut)
                    {
                        *pmask = (*func)(r * rc_inv);
                        assert(*pmask <= 1);
                        assert(*pmask >= 0);
                        ncount++;
                        assert(icount);
                    }
                    else
                    {
                        *pmask = 0;
                    }
                }
#ifdef DEBUG
                lcount++;
                assert(ix * incx + iy * incy + iz < size);
#endif

                xc[2] += hgrid[2];
                pmask++;
            }
            xc[1] += hgrid[1];
        }
        xc[0] += hgrid[0];
    }

    // icount needs to be at least 1 so that the next check for overlap is
    // consistent (e.g. in restart) even if no actual point was found in
    // sphere!!! Otherwise restart may attribute data to wrong orbital when
    // changing radius
    if (level == 0) icount = std::max(1, icount);

#if 0
    if( icount==0 && level==0){
        const double r = sub_cell_[level]->distance(center, iloc, ct.bcPoisson);
        (*MPIdata::sout)<<"delta="<<delta<<", r="<<r<<", Centers :"<<center<<"---";
        sub_cell_[level]->printCenter(iloc,(*MPIdata::sout));
    }
#endif

    assert((pmask - 1) == &mask[size - 1]);

    assign(iloc, level, ncount, icount, mask);

#pragma omp master
    init_tm_.stop();

    return icount;
}

void GridMask::plot(const unsigned short level, const int st)
{
    char filename[20], extension[20];
    sprintf(filename, "%s", "mask_");
    sprintf(extension, "%d", st);
    strcat(filename, extension);
    sprintf(extension, "%s", "_");
    strcat(filename, extension);
    sprintf(extension, "%d", level);
    strcat(filename, extension);

    std::ofstream ofile;
    ofile.open(&filename[0]);
    const int lnumpt = loc_numpt_[level];

    for (unsigned short iloc = 0; iloc < subdivx_; iloc++)
    {
        if (mask_not_zero_[level][iloc] == 0)
        {
            for (int i = 0; i < lnumpt; i++)
                ofile << 0 << std::endl;
        }
        else if (mask_not_zero_[level][iloc] == 1)
        {
            for (int i = 0; i < lnumpt; i++)
                ofile << 1 << std::endl;
        }
        else if (mask_not_zero_[level][iloc] == 2)
        {
            lmasktype* pmask = &lmask_[level][iloc][0];
            for (int i = 0; i < lnumpt; i++)
                ofile << (int)pmask[i] << std::endl;
        }
    }
    ofile.close();
}

template <typename T>
void GridMask::apply(pb::GridFunc<T>& gu, const unsigned short level)
{
    for (unsigned short iloc = 0; iloc < subdivx_; iloc++)
        apply(gu, level, iloc);
}

void GridMask::assign(const unsigned short iloc, const unsigned short level,
    const int num, const int xnum, const std::vector<lmasktype>& val)
{
    assert(xnum >= num);
    assert(iloc < subdivx_);
    assert(level < (int)mask_not_zero_.size());
    assert((int)val.size() >= num);
    assert((int)val.size() >= xnum);
    assert((int)lmask_[level].size() > iloc);

    // Early return if all elements are 0 and stay 0
    if (num == 0 && xnum == 0)
    {
        if (mask_not_zero_[level][iloc] == 2) lmask_[level][iloc].resize(0);
        mask_not_zero_[level][iloc] = -1;
        return;
    }

    // Early return if all elements are 0
    if (num == 0)
    {
        if (mask_not_zero_[level][iloc] == 2) lmask_[level][iloc].resize(0);
        mask_not_zero_[level][iloc] = 0;
        return;
    }
    mask_not_zero_[level][iloc] = 2;

    lmask_[level][iloc].resize(val.size());
    lmask_[level][iloc] = val;
    assert(lmask_[level][iloc].size() == val.size());
}

void GridMask::assign1(const unsigned short iloc, const unsigned short level)
{
    assert(iloc < subdivx_);
    assert(level < (unsigned short)mask_not_zero_.size());
    assert((unsigned short)lmask_[level].size() > iloc);

    if (mask_not_zero_[level][iloc] == 2) lmask_[level][iloc].resize(0);
    mask_not_zero_[level][iloc] = 1;

    return;
}

bool GridMask::overlap_on_pe(
    const GridMask& gm, const unsigned short level) const
{
    assert(level < mask_not_zero_.size());
    assert(level < gm.mask_not_zero_.size());
    assert(subdivx_ <= mask_not_zero_[level].size());
    assert(subdivx_ <= gm.mask_not_zero_[level].size());

    bool overlap_12 = false;

    // sum up the contributions of subdomains
    for (unsigned short iloc = 0; iloc < subdivx_; iloc++)
    {
        assert(mask_not_zero_[level][iloc] < 3);
        assert(gm.mask_not_zero_[level][iloc] < 3);
        assert(mask_not_zero_[level][iloc] >= -1);
        assert(gm.mask_not_zero_[level][iloc] >= -1);
        overlap_12 = (overlap_12
                      || ((mask_not_zero_[level][iloc] + 1)
                          && (gm.mask_not_zero_[level][iloc] + 1)));
    }

    return overlap_12;
}

template <typename T>
void GridMask::multiplyByMask(
    T* u, const unsigned short level, const unsigned short iloc) const
{
    const int lnumpt = loc_numpt_[level];
    const std::vector<lmasktype>& maskptr(lmask_[level][iloc]);
    T* const pu = u + offset(level, iloc);
    for (int idx = 0; idx < lnumpt; idx++)
    {
        assert(maskptr[idx] >= 0 && maskptr[idx] <= 1);
        pu[idx] *= (T)maskptr[idx];
    }
}

template <typename T>
void GridMask::cutWithMask(
    T* u, const unsigned short level, const unsigned short iloc) const
{
    const int lnumpt = loc_numpt_[level];
    const std::vector<lmasktype>& maskptr(lmask_[level][iloc]);
    T* const pu = u + offset(level, iloc);
    for (int idx = 0; idx < lnumpt; idx++)
    {
        assert(maskptr[idx] >= 0 && maskptr[idx] <= 1);
        pu[idx] = limitAbsValue(pu[idx], maskptr[idx]);
    }
}

// explicit instantiations
template void GridMask::apply(
    pb::GridFunc<ORBDTYPE>& gu, const unsigned short level);

template void GridMask::multiplyByMask(
    float* u, const unsigned short level, const unsigned short iloc) const;
template void GridMask::multiplyByMask(
    double* u, const unsigned short level, const unsigned short iloc) const;

template void GridMask::cutWithMask(
    float* u, const unsigned short level, const unsigned short iloc) const;
template void GridMask::cutWithMask(
    double* u, const unsigned short level, const unsigned short iloc) const;
