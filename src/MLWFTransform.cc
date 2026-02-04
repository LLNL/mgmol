// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// #define DEBUG 1
#include "MGmol_blas1.h"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Control.h"
#include "MGmol_MPI.h"
#include "MLWFTransform.h"

#include <mpi.h>

using namespace std;

double jade(std::vector<dist_matrix::DistMatrix<double>*>& r_, const int n_loc,
    const int m_loc, const int offset_, const int maxsweep, const double tol,
    MPI_Comm comm, std::vector<double>& mat_, const bool print_flag);

MLWFTransform::MLWFTransform(
    const int nst, const Vector3D& origin, const Vector3D& ll)
    : OrbitalsTransform(nst, origin, ll)
{
    const int npcol = bcr_->npcol();
    nstcol_         = bsize_ * npcol;

    if (nst > 0)
    {
        std::vector<std::string> names(2 * NDIM);
        names[0] = "matrixR_0";
        names[1] = "matrixR_1";
        names[2] = "matrixR_2";
        names[3] = "matrixR_3";
#if NDIM > 2
        names[4] = "matrixR_4";
        names[5] = "matrixR_5";
#endif
        for (int k = 0; k < 2 * NDIM; k++)
        {
            r_[k] = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                names[k], *bcr_, nst_, nstcol_, nst_, bsize_);
        }
    }
}

double MLWFTransform::spread2(int i, int j) const
{
    assert(i >= 0 && i < nst_);
    assert(j >= 0 && j < NDIM);
    assert(cell_[j] > 0.);

    int ii = -1;
    double cs[2];

    // integral over domain of position operator^2 (sin^2+cos^2) * 0.5 * M_1_PI
    const double lby2pi = 0.5 * M_1_PI * cell_[j];
    if (offset_ >= 0 && i >= offset_ && i < offset_ + bsize_)
    {
        ii = nst_ * (i - offset_) + i;
        assert(ii >= 0);
        cs[0] = r_[2 * j]->val(ii);
        cs[1] = r_[2 * j + 1]->val(ii);
    }
    else
    {
        cs[0] = 0.;
        cs[1] = 0.;
    }

    if (offset_ > 0 && ii != -1)
    {
        MPI_Send(&cs[0], 2, MPI_DOUBLE, 0, 0, comm_);
    }
    if (offset_ == 0 && ii == -1)
    {
        MPI_Status status;
        int src = i / bsize_;
        MPI_Recv(&cs[0], 2, MPI_DOUBLE, src, 0, comm_, &status);
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.bcast(cs, 2);

    return lby2pi * lby2pi - (cs[0] * cs[0] + cs[1] * cs[1]);
}

////////////////////////////////////////////////////////////////////////////////

void MLWFTransform::setia(vector<int>& iiu) { iu_ = iiu; }

////////////////////////////////////////////////////////////////////////////////

void MLWFTransform::compute_transform(const int maxsweep, const double tol)
{
    if (bcr_->onpe0())
        (*MPIdata::sout) << " MLWFTransform for nst=" << nst_ << endl;

    // sine and cosine operators only: compute MLWF transform

    double delta = jade(
        r_, bsize_, nst_, offset_, maxsweep, tol, comm_, mat_, bcr_->onpe0());

    if (bcr_->onpe0())
        if (delta > tol)
            (*MPIdata::sout)
                << " MLWFTransform: decrease was " << setprecision(12) << delta
                << " after " << maxsweep << " iterations" << endl;
}

void MLWFTransform::printTransform()
{
    Control& ct     = *(Control::instance());
    const int numst = ct.numst;

    if (onpe0)
    {
        filebuf fb;
        fb.open("transformMLWF.dat", ios::out);
        ostream os(&fb);
        os << numst << endl << endl;

        for (int j = 0; j < numst; j++)
        {
            for (int i = 0; i < numst; i++)
                os << mat_[j * numst + i] << endl;
            os << endl;
        }
        fb.close();
    }
}
// #endif
