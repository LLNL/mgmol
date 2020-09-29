// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "OrbitalsTransform.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include <iomanip>

OrbitalsTransform::OrbitalsTransform(
    const int nst, const Vector3D& origin, const Vector3D& ll)
    : nst_(nst), cell_(ll), origin_(origin)
{
    assert(nst >= 0);
    assert(nst < 100000);
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();
    bcr_  = new dist_matrix::BlacsContext(comm, 'r', 1 + (nst - 1) / 4);
    comm_ = bcr_->comm_active();
    const int npcol = bcr_->npcol();

    const int ndir = 6;
    r_.resize(ndir, nullptr);

    bsize_ = nst_ / npcol;
    if (bsize_ * bcr_->npcol() < nst_) bsize_++;
    if (bsize_ % 2 != 0) bsize_++;

    if (onpe0)
        (*MPIdata::sout) << " OrbitalsTransform using " << bcr_->npcol()
                         << " parallel tasks,"
                         << " with " << bsize_ << " columns per task"
                         << std::endl;

    if (bcr_->active())
    {
        const int mypecol = bcr_->mycol();
        offset_           = mypecol * bsize_;
        lnst_             = std::min(bsize_, nst_ - offset_);
        if (lnst_ < 0) lnst_ = 0;
        iu_.resize(nst_ * bsize_, 1);
    }
    else
    {
        lnst_   = 0;
        offset_ = -1;
    }
    mat_.resize(nst_ * nst_);

    assert(lnst_ <= nst);
    assert(lnst_ >= 0);
}

OrbitalsTransform::~OrbitalsTransform()
{
    for (int k = 0; k < 2 * NDIM; k++)
    {
        delete r_[k];
    }
    delete bcr_;
}

////////////////////////////////////////////////////////////////////////////////

void OrbitalsTransform::printTransformationMatrix() const
{
    if (onpe0)
    {

        (*MPIdata::sout) << "Transform matrix: " << std::endl;
        (*MPIdata::sout) << std::setprecision(5);
        for (int i = 0; i < nst_; ++i)
        {
            for (int j = 0; j < nst_; ++j)
                (*MPIdata::sout) << mat_[j * nst_ + i] << "  ";
            (*MPIdata::sout) << std::endl;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

Vector3D OrbitalsTransform::center(void) const
{
    Vector3D sum(0.0, 0.0, 0.0);
    for (int i = 0; i < nst_; i++)
        sum += center(i);
    return sum;
}

////////////////////////////////////////////////////////////////////////////////
// return MLWC; reference=origin of domain at (0,0,0)

Vector3D OrbitalsTransform::center(int i) const
{
    assert(i >= 0 && i < nst_);

    Vector3D v;
    int ii;
    double cx, sx, cy, sy;
#if (NDIM == 3)
    double cz, sz;
#endif

    if (i >= offset_ && i < offset_ + bsize_ && offset_ >= 0)
    {
        ii = nst_ * (i - offset_) + i;
        cx = r_[0]->val(ii);
        sx = r_[1]->val(ii);
        cy = r_[2]->val(ii);
        sy = r_[3]->val(ii);
#if (NDIM == 3)
        cz = r_[4]->val(ii);
        sz = r_[5]->val(ii);
#endif
    }
    else
        ii = -1;

    if (ii != -1)
    {
        const double itwopi = 1.0 / (2.0 * M_PI);
        // get center between -0.5*cell and +0.5*cell
        v[0] = cell_[0] * itwopi * atan2(sx, cx);
        v[1] = cell_[1] * itwopi * atan2(sy, cy);
#if (NDIM == 3)
        v[2] = cell_[2] * itwopi * atan2(sz, cz);
#endif
        // get center between 0. and cell
        for (int i = 0; i < NDIM; i++)
            if (v[i] < 0.) v[i] += cell_[i];
        // get center in cell
        for (int i = 0; i < NDIM; i++)
            v[i] += origin_[i];
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Status status;
    if (offset_ > 0 && ii != -1) MPI_Send(&v[0], NDIM, MPI_DOUBLE, 0, 0, comm_);
    if (offset_ == 0 && ii == -1)
    {
        int src = i / bsize_;
        MPI_Recv(&v[0], NDIM, MPI_DOUBLE, src, 0, comm_, &status);
    }
    mmpi.bcast(&v[0], NDIM);
    return v;
}

////////////////////////////////////////////////////////////////////////////////
double OrbitalsTransform::spread2(int i) const
{
    assert((i >= 0) && (i < nst_));
    return spread2(i, 0) + spread2(i, 1) + spread2(i, 2);
}

////////////////////////////////////////////////////////////////////////////////
double OrbitalsTransform::spread(int i) const
{
    assert(spread2(i) > 0.);
    return sqrt(spread2(i));
}

////////////////////////////////////////////////////////////////////////////////
double OrbitalsTransform::spread2(void) const
{
    double sum = 0.0;
    for (int i = 0; i < nst_; i++)
        sum += spread2(i);
    return sum;
}

////////////////////////////////////////////////////////////////////////////////
double OrbitalsTransform::spread(void) const { return std::sqrt(spread2()); }

////////////////////////////////////////////////////////////////////////////////
double OrbitalsTransform::volume() const
{
    double vol = 0.;
    for (int i = 0; i < nst_; i++)
    {
        const double spreadi = spread(i);
        vol += spreadi * spreadi * spreadi;
    }
    return 4. * vol * M_PI / 3.;
}

////////////////////////////////////////////////////////////////////////////////
void OrbitalsTransform::printCentersAndSpreads(std::ostream& os)
{
    if (bcr_->onpe0() && nst_ > 0)
    {
        os << "------------------------------------------------------"
           << std::endl;
        os << " Orbitals centers and spreads " << std::endl << std::endl;
    }
    double max = 0.;
    double min = 1.e12;

    for (int i = 0; i < nst_; i++)
    {
        Vector3D centeri(center(i));
        double spreadi = spread(i);
        if (bcr_->onpe0())
        {
            os << "&& " << std::setw(4) << i + 1 << "   ";
            os.setf(std::ios::fixed, std::ios::floatfield);
            os.setf(std::ios::right, std::ios::adjustfield);
            os << std::setw(12) << std::setprecision(3) << centeri[0] << " "
               << std::setw(12) << std::setprecision(3) << centeri[1] << " "
               << std::setw(12) << std::setprecision(3) << centeri[2]
               << "         ";
            os << std::setw(12) << std::setprecision(3) << spreadi << std::endl;
        }
        if (spreadi > max) max = spreadi;
        if (spreadi < min) min = spreadi;
    }

    if (bcr_->onpe0())
        os << "Min. spread = " << min << ", Max. spread = " << max << std::endl
           << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
void OrbitalsTransform::distributeColumnsR(
    std::vector<std::vector<DISTMATDTYPE>>& vmm)
{
    assert(vmm.size() == 6);

    // mycol =0 if mat not distributed (duplicated)
    const int mycol = r_[0]->mycol();

    if (mycol < 0) return;

    // calculate number of elements to assign
    int count = bsize_ * nst_; // size of matrix block
    if (static_cast<int>(vmm[0].size()) < (mycol + 1) * bsize_ * nst_)
        count = vmm[0].size() - mycol * bsize_ * nst_;

    if (count > 0)
    {
        for (int i = 0; i < 6; i++)
        {
            assert(static_cast<int>(vmm[i].size()) > mycol * bsize_ * nst_);
            r_[i]->assign(&vmm[i][mycol * bsize_ * nst_], count);
        }
    }
}
