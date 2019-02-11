// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "global.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "ColoredRegions.h"
#include "Control.h"
#include "DistMatrix.h"
#include "FunctionsPacking.h"
#include "GridFunc.h"
#include "GridMask.h"
#include "HDFrestart.h"
#include "Laph2.h"
#include "Laph4M.h"
#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "MPIdata.h"
#include "Masks4Orbitals.h"
#include "MasksSet.h"
#include "Mesh.h"
#include "Potentials.h"
#include "Preconditioning.h"
#include "ProjectedMatrices.h"
#include "ReplicatedWorkSpace.h"
#include "SparseDistMatrix.h"
#include "SquareLocalMatrices.h"
#include "SubCell.h"
#include "SubMatrices.h"
#include "VariableSizeMatrix.h"
#include "hdf_tools.h"
#include "lapack_c.h"

#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

#define ORBITAL_OCCUPATION 2.
string getDatasetName(const string& name, const int color);

short LocGridOrbitals::subdivx_          = 0;
int LocGridOrbitals::lda_                = 0;
int LocGridOrbitals::numpt_              = 0;
int LocGridOrbitals::loc_numpt_          = 0;
short LocGridOrbitals::bc_[3]            = { 1, 1, 1 };
PtrFunc LocGridOrbitals::dotProduct_     = &LocGridOrbitals::dotProductDiagonal;
int LocGridOrbitals::data_wghosts_index_ = -1;

Timer LocGridOrbitals::get_dm_tm_("LocGridOrbitals::get_dm");
Timer LocGridOrbitals::matB_tm_("LocGridOrbitals::matB");
Timer LocGridOrbitals::invBmat_tm_("LocGridOrbitals::invBmat");
Timer LocGridOrbitals::overlap_tm_("LocGridOrbitals::overlap");
Timer LocGridOrbitals::dot_product_tm_("LocGridOrbitals::dot_product");
Timer LocGridOrbitals::addDot_tm_("LocGridOrbitals::addDot");
Timer LocGridOrbitals::mask_tm_("LocGridOrbitals::mask");
Timer LocGridOrbitals::prod_matrix_tm_("LocGridOrbitals::prod_matrix");
Timer LocGridOrbitals::assign_tm_("LocGridOrbitals::assign");
Timer LocGridOrbitals::normalize_tm_("LocGridOrbitals::normalize");
Timer LocGridOrbitals::axpy_tm_("LocGridOrbitals::axpy");

LocGridOrbitals::LocGridOrbitals(std::string name,
    const pb::Grid& my_grid, const short subdivx,
    const int numst, const short bc[3],
    ProjectedMatricesInterface* proj_matrices, LocalizationRegions* lrs,
    MasksSet* masks, MasksSet* corrmasks, ClusterOrbitals* local_cluster,
    const bool setup_flag)
    : name_(name),
      grid_(my_grid),
      proj_matrices_(proj_matrices),
      block_vector_(my_grid, subdivx, bc),
      lrs_(lrs),
      local_cluster_(local_cluster)
{
    // preconditions
    assert(subdivx > 0);
    assert(proj_matrices != 0);
    assert(lrs != 0);

    for (short i = 0; i < 3; i++)
        assert(bc[i] == 0 || bc[i] == 1);
    assert(grid_.size() > 0);

    subdivx_   = subdivx;
    numst_     = numst;
    numpt_     = grid_.size();
    lda_       = block_vector_.getld();
    loc_numpt_ = numpt_ / subdivx_;

    chromatic_number_ = 0;

    assert(numst_ >= 0);

    gidToStorage_ = 0;

    overlapping_gids_.clear();

    if (masks && corrmasks)
        masks4orbitals_.reset(
            new Masks4Orbitals(masks, corrmasks, lrs->getOverlapGids()));

    if (setup_flag) setup(lrs);
}

LocGridOrbitals::~LocGridOrbitals()
{
    assert(proj_matrices_ != 0);
    assert(pack_);
    assert(gidToStorage_ != 0);

    // delete gidToStorage here. This is OK since it is not a shared data
    // else there would be a memory leak.
    delete gidToStorage_;
    gidToStorage_ = 0;
}

LocGridOrbitals::LocGridOrbitals(const std::string name,
    const LocGridOrbitals& A, const bool copy_data)
    : Orbitals(A, copy_data),
      name_(name),
      grid_(A.grid_),
      proj_matrices_(A.proj_matrices_),
      block_vector_(A.block_vector_, copy_data),
      masks4orbitals_(A.masks4orbitals_),
      lrs_(A.lrs_),
      local_cluster_(A.local_cluster_)
{
    // if(onpe0)cout<<"call LocGridOrbitals(const LocGridOrbitals &A, const bool
    // copy_data)"<<endl;

    assert(A.chromatic_number_ >= 0);
    assert(A.proj_matrices_ != 0);
    assert(A.lrs_ != 0);

    copySharedData(A);

    gidToStorage_ = 0;

    setGids2Storage();
}

LocGridOrbitals::LocGridOrbitals(const std::string name,
    const LocGridOrbitals& A,
    ProjectedMatricesInterface* proj_matrices, MasksSet* masks,
    MasksSet* corrmasks, const bool copy_data)
    : Orbitals(A, copy_data),
      name_(name),
      grid_(A.grid_),
      proj_matrices_(proj_matrices),
      block_vector_(A.block_vector_, copy_data),
      lrs_(A.lrs_),
      local_cluster_(A.local_cluster_)
{
    assert(A.chromatic_number_ >= 0);
    assert(proj_matrices != 0);
    assert(masks != 0);
    assert(lrs_ != 0);

    copySharedData(A);

    gidToStorage_ = 0;

    setGids2Storage();

    masks4orbitals_.reset(
        new Masks4Orbitals(masks, corrmasks, A.all_overlapping_gids_));

    // setup new projected_matrices object
    Control& ct = *(Control::instance());
    proj_matrices_->setup(ct.occ_width, ct.getNel(), overlapping_gids_);
}

void LocGridOrbitals::copySharedData(const LocGridOrbitals& A)
{
    // if(onpe0)cout<<"call LocGridOrbitals::copySharedData(const
    // LocGridOrbitals &A)"<<endl;

    assert(A.gidToStorage_ != 0);
    assert(A.pack_);

    numst_ = A.numst_;

    lrs_iterative_index_ = A.lrs_iterative_index_;

    chromatic_number_ = A.chromatic_number_;

    pack_ = A.pack_;

    overlapping_gids_     = A.overlapping_gids_;
    all_overlapping_gids_ = A.all_overlapping_gids_;
    local_gids_           = A.local_gids_;

    distributor_diagdotprod_ = A.distributor_diagdotprod_;
    distributor_normalize_   = A.distributor_normalize_;
}

void LocGridOrbitals::copyDataFrom(const LocGridOrbitals& src)
{
    assert(proj_matrices_ != 0);

    block_vector_.copyDataFrom(src.block_vector_);

    setIterativeIndex(src);
}

void LocGridOrbitals::setDotProduct(const short dot_type)
{
    if (dot_type == 0)
        dotProduct_ = &LocGridOrbitals::dotProductDiagonal;
    else if (dot_type == 1)
        dotProduct_ = &LocGridOrbitals::dotProductWithInvS;
    else if (dot_type == 2)
        dotProduct_ = &LocGridOrbitals::dotProductWithDM;
    else if (dot_type == 3)
        dotProduct_ = &LocGridOrbitals::dotProductSimple;
}

void LocGridOrbitals::setGids2Storage()
{
    assert(chromatic_number_ >= 0);
    assert(subdivx_ > 0);

    if (gidToStorage_ != 0)
        gidToStorage_->clear();
    else
        gidToStorage_ = new vector<map<int, ORBDTYPE*>>();
    gidToStorage_->resize(subdivx_);
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        map<int, ORBDTYPE*>& gid2st((*gidToStorage_)[iloc]);
        for (int color = 0; color < chromatic_number_; color++)
        {
            const int gid = overlapping_gids_[iloc][color];
            if (gid != -1)
            {
                gid2st.insert(pair<int, ORBDTYPE*>(gid, getPsi(color, iloc)));
            }
        }
    }
}

// return pointer to const data
const ORBDTYPE* LocGridOrbitals::getGidStorage(
    const int gid, const short iloc) const
{
    assert(numst_ >= 0);
    assert(iloc < subdivx_);
    assert(iloc >= 0);
    assert(gid < numst_);
    assert(iloc < (short)gidToStorage_->size());

    map<int, ORBDTYPE*>::const_iterator p = (*gidToStorage_)[iloc].find(gid);
    if (p != (*gidToStorage_)[iloc].end())
        return p->second;
    else
        return 0;
}

void LocGridOrbitals::setup(
    MasksSet* masks, MasksSet* corrmasks, LocalizationRegions* lrs)
{
    assert(masks != 0);
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp("LocGridOrbitals::setup(MasksSet*, MasksSet*)...",
            (*MPIdata::sout));

    masks4orbitals_.reset(
        new Masks4Orbitals(masks, corrmasks, lrs->getOverlapGids()));

    setup(lrs);
}

void LocGridOrbitals::setup(LocalizationRegions* lrs)
{
    Control& ct = *(Control::instance());

    // preconditions
    assert(lrs != 0);
    assert(proj_matrices_ != 0);

    if (ct.verbose > 0)
        printWithTimeStamp("LocGridOrbitals::setup()...", (*MPIdata::sout));

    lrs_iterative_index_ = lrs->getIterativeIndex();

    overlapping_gids_.clear();

    chromatic_number_ = packStates(lrs);

    computeGlobalIndexes(*lrs);

    bool skinny_stencil = !ct.Mehrstellen();

    block_vector_.initialize(overlapping_gids_, skinny_stencil);

    setGids2Storage();

    proj_matrices_->setup(ct.occ_width, ct.getNel(), overlapping_gids_);

    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };

    double maxr = lrs->max_radii();
    distributor_diagdotprod_.reset(
        new DataDistribution("dot", maxr, myPEenv, domain));
    distributor_normalize_.reset(
        new DataDistribution("norm", 2. * maxr, myPEenv, domain));

    if (ct.verbose > 0)
        printWithTimeStamp(
            "LocGridOrbitals::setup() done...", (*MPIdata::sout));
}

void LocGridOrbitals::reset(
    MasksSet* masks, MasksSet* corrmasks, LocalizationRegions* lrs)
{
    // free some old data
    block_vector_.clear();
    setIterativeIndex(-10);

    // reset
    setup(masks, corrmasks, lrs);
}

void LocGridOrbitals::assign(const LocGridOrbitals& orbitals)
{
    assign_tm_.start();

#ifdef DEBUG
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.barrier();
#endif
    assert(chromatic_number_ >= 0);
    assert(proj_matrices_ != 0);

    setIterativeIndex(orbitals);

    // if( onpe0 && ct.verbose>2 )
    //    (*MPIdata::sout)<<"LocGridOrbitals::Assign orbitals up to
    //    "<<chromatic_number_<<endl;

    if (pack_ == orbitals.pack_)
    {

        block_vector_.copyDataFrom(orbitals.block_vector_);
    }
    else
    {
        Control& ct = *(Control::instance());
        if (onpe0 && ct.verbose > 2)
            (*MPIdata::sout)
                << "LocGridOrbitals::Assign orbitals to different LR" << endl;
        for (int color = 0; color < chromatic_number_; color++)
        {
            // assign state
            for (short iloc = 0; iloc < subdivx_; iloc++)
            {
                const int gid = overlapping_gids_[iloc][color];
                if (gid != -1)
                {
                    // find storage location in orbitals
                    const ORBDTYPE* const val
                        = orbitals.getGidStorage(gid, iloc);
                    // copy into new psi_
                    if (val != 0)
                    {
                        block_vector_.assignLocal(color, iloc, val);
                    }
                }
            }
        }
    }

    assign_tm_.stop();
}

void LocGridOrbitals::axpy(const double alpha, const LocGridOrbitals& orbitals)
{
    axpy_tm_.start();

    assert(pack_ != 0);
    assert(orbitals.pack_ != 0);
    assert(overlapping_gids_.size() > 0);

    //    int ione=1;

    if (pack_ == orbitals.pack_)
    {
        block_vector_.axpy(alpha, orbitals.block_vector_);
    }
    else
    {
        for (int color = 0; color < chromatic_number_; color++)
        {
            // assign state
            for (short iloc = 0; iloc < subdivx_; iloc++)
            {
                const int gid = overlapping_gids_[iloc][color];
                assert(gid < numst_);
                if (gid != -1)
                {
                    // find orbital storage in orbitals
                    const ORBDTYPE* const val
                        = orbitals.getGidStorage(gid, iloc);
                    // copy into new psi_
                    if (val != 0)
                    {
                        MPaxpy(loc_numpt_, alpha, val, getPsi(color, iloc));
                    }
                }
            }
        }
    }
    incrementIterativeIndex();

    axpy_tm_.stop();
}

short LocGridOrbitals::checkOverlap(
    const int st1, const int st2, const short level)
{
    assert(masks4orbitals_);

    return masks4orbitals_->checkOverlap(st1, st2, level);
}

void LocGridOrbitals::applyMask(const bool first_time)
{
    assert(chromatic_number_ >= 0);
    assert(subdivx_ > 0);
    assert(loc_numpt_ > 0);

    mask_tm_.start();

    const short ln = 0;
    for (int color = 0; color < chromatic_number_; color++)
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const int gid = overlapping_gids_[iloc][color];
            if (gid != -1)
            {
                (masks4orbitals_->getMask(gid))
                    .apply(psi(color), ln, iloc, first_time);
            }
            else
                block_vector_.set_zero(color, iloc);
        }

    incrementIterativeIndex();

    mask_tm_.stop();
}

void LocGridOrbitals::applyCorrMask(const bool first_time)
{
    mask_tm_.start();

    for (int color = 0; color < chromatic_number_; color++)
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const int gid = overlapping_gids_[iloc][color];
            if (gid != -1)
            {
                (masks4orbitals_->getCorrMask(gid))
                    .apply(psi(color), 0, iloc, first_time);
            }
            else
                block_vector_.set_zero(color, iloc);
        }
    incrementIterativeIndex();

    mask_tm_.stop();
}

void LocGridOrbitals::app_mask(
    const int color, ORBDTYPE* u, const short level) const
{
    mask_tm_.start();
    assert(masks4orbitals_);

    assert(color < chromatic_number_);
    int lnumpt = (loc_numpt_ >> level);
    lnumpt     = (lnumpt >> level);
    lnumpt     = (lnumpt >> level);

    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        int gid = overlapping_gids_[iloc][color];
        if (gid != -1)
        {
            (masks4orbitals_->getMask(gid)).apply(u, level, iloc);
        }
        else
            memset(u + iloc * lnumpt, 0, lnumpt * sizeof(ORBDTYPE));
    }
    mask_tm_.stop();
}

void LocGridOrbitals::app_mask(
    const int color, pb::GridFunc<ORBDTYPE>& gu, const short level) const
{
    mask_tm_.start();

    assert(color < chromatic_number_);
    assert(masks4orbitals_);

    const short shift = gu.grid().ghost_pt();
    const int dim0    = gu.grid().dim(0) / subdivx_;
    const int incx    = gu.grid().inc(0);

    const int lnumpt = dim0 * incx;

    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        const int gid = overlapping_gids_[iloc][color];

        if (gid != -1)
        {
            (masks4orbitals_->getMask(gid)).apply(gu, level, iloc);
        }
        else
        {
            int offset = (shift + dim0 * iloc) * incx;
            assert(offset + lnumpt < gu.grid().sizeg());
            ORBDTYPE* pu = gu.uu() + offset;
            memset(pu, 0, lnumpt * sizeof(ORBDTYPE));
        }
    }
    mask_tm_.stop();
}

void LocGridOrbitals::init2zero()
{
    for (int icolor = 0; icolor < chromatic_number_; icolor++)
    {
        ORBDTYPE* ipsi = psi(icolor);
        memset(ipsi, 0, numpt_ * sizeof(ORBDTYPE));
    }
}

void LocGridOrbitals::initGauss(const double rc, const LocalizationRegions& lrs)
{
    assert(chromatic_number_ >= 0);
    assert(subdivx_ > 0);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Control& ct     = *(Control::instance());
    if (mmpi.instancePE0() && ct.verbose > 2)
        (*MPIdata::sout) << "Initial orbitals: Gaussians of width " << rc
                         << endl;

    const double invrc2 = 1. / (rc * rc);

    const double start0 = grid_.start(0);
    const double start1 = grid_.start(1);
    const double start2 = grid_.start(2);

    const int dim0 = grid_.dim(0) / subdivx_;
    const int dim1 = grid_.dim(1);
    const int dim2 = grid_.dim(2);

    const int incx = dim1 * dim2;
    const int incy = dim2;

    const double hgrid[3] = { grid_.hgrid(0), grid_.hgrid(2), grid_.hgrid(2) };

    Vector3D ll;
    for (short i = 0; i < 3; i++)
        ll[i] = grid_.ll(i);

    const double rmax = 6. * rc;
    for (int icolor = 0; icolor < chromatic_number_; icolor++)
    {
        ORBDTYPE* ipsi = psi(icolor);
        memset(ipsi, 0, numpt_ * sizeof(ORBDTYPE));

        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const int gid = overlapping_gids_[iloc][icolor];
            if (gid > -1)
            {
                const Vector3D& center(lrs.getCenter(gid));
                Vector3D xc;

                xc[0] = start0 + iloc * dim0 * hgrid[0];
                for (int ix = iloc * dim0; ix < (iloc + 1) * dim0; ix++)
                {
                    xc[1] = start1;

                    for (int iy = 0; iy < dim1; iy++)
                    {
                        xc[2] = start2;
                        for (int iz = 0; iz < dim2; iz++)
                        {
                            const double r = xc.minimage(center, ll, bc_);
                            if (r < rmax)
                                ipsi[ix * incx + iy * incy + iz]
                                    = (ORBDTYPE)exp(-r * r * invrc2);
                            else
                                ipsi[ix * incx + iy * incy + iz] = 0.;

                            xc[2] += hgrid[2];
                        }
                        xc[1] += hgrid[1];
                    }
                    xc[0] += hgrid[0];
                }
            }
        }
    }
    resetIterativeIndex();
}

void LocGridOrbitals::initFourier()
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << "Initial orbitals: Fourier " << endl;

    const double start0 = grid_.start(0);
    const double start1 = grid_.start(1);
    const double start2 = grid_.start(2);

    const int dim0 = grid_.dim(0) / subdivx_;
    const int dim1 = grid_.dim(1);
    const int dim2 = grid_.dim(2);

    const int incx = dim1 * dim2;
    const int incy = dim2;

    const double hgrid[3] = { grid_.hgrid(0), grid_.hgrid(1), grid_.hgrid(2) };

    Vector3D ll;
    for (short i = 0; i < 3; i++)
        ll[i] = grid_.ll(i);

    const double dk[3]
        = { 2. * M_PI / ll[0], 2. * M_PI / ll[1], 2. * M_PI / ll[2] };

    const int cbrtncolors = (int)ceil(cbrt(chromatic_number_));

    for (int icolor = 0; icolor < chromatic_number_; icolor++)
    {
        int k0 = icolor / (cbrtncolors * cbrtncolors);
        int k1 = (icolor - k0 * cbrtncolors * cbrtncolors) / cbrtncolors;
        int k2 = icolor % cbrtncolors;
        // if( onpe0 )(*MPIdata::sout)<<" k=("<<k0<<","<<k1<<","<<k2<<")"<<endl;

        const double kk[3]
            = { dk[0] * (double)k0, dk[1] * (double)k1, dk[2] * (double)k2 };

        ORBDTYPE* ipsi = psi(icolor);
        memset(ipsi, 0, numpt_ * sizeof(ORBDTYPE));

        for (short iloc = 0; iloc < subdivx_; iloc++)
        {

            const int gid = overlapping_gids_[iloc][icolor];
            if (gid > -1)
            {
                double x = start0 + iloc * dim0 * hgrid[0];
                for (int ix = iloc * dim0; ix < (iloc + 1) * dim0; ix++)
                {
                    double y = start1;

                    for (int iy = 0; iy < dim1; iy++)
                    {
                        double z = start2;
                        for (int iz = 0; iz < dim2; iz++)
                        {
                            ipsi[ix * incx + iy * incy + iz]
                                = (ORBDTYPE)(cos(kk[0] * x) * cos(kk[1] * y)
                                             * cos(kk[2] * z));

                            z += hgrid[2];
                        }
                        y += hgrid[1];
                    }
                    x += hgrid[0];
                }
            }
        }
    }
    resetIterativeIndex();
}

int LocGridOrbitals::packStates(LocalizationRegions* lrs)
{
    assert(lrs != 0);

    Control& ct = *(Control::instance());

    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << " PACK STATES " << 0 << " to " << numst_ << endl;

    // compute overlap for all the orbitals on all PEs
    const int dim = ct.globalColoring() ? numst_ : lrs->getNumOverlapGids();

    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << " PACK " << dim << " STATES" << endl;

    const bool global = ct.globalColoring();

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    pack_.reset(new FunctionsPacking(lrs, global, mmpi.commSameSpin()));

    assert(pack_->chromatic_number() < 100000);

    return pack_->chromatic_number();
}

void LocGridOrbitals::precond_smooth(ORBDTYPE* rhs, const int ld,
    const int ifirst, const int nvect, const int npower, const double alpha)
{
    assert(ld >= grid_.size());

    pb::Laph2<ORBDTYPE> myoper(grid_);
    pb::GridFunc<ORBDTYPE> gf_w(grid_, bc_[0], bc_[1], bc_[2]);

    for (int i = 0; i < nvect; i++)
    {
        ORBDTYPE* rhsi = rhs + ld * i;

        pb::GridFunc<ORBDTYPE> gf_sd(rhsi, grid_, bc_[0], bc_[1], bc_[2]);

        myoper.smooth(gf_sd, gf_w, alpha);

        app_mask(ifirst + i, gf_w, 0);

        gf_w.init_vect(rhsi, 'd');
    }
}

void LocGridOrbitals::multiply_by_matrix(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& dmatrix,
    ORBDTYPE* const product, const int ldp)
{
#if 0
    (*MPIdata::sout)<<"self multiply_by_matrix"<<endl;
#endif

    ReplicatedWorkSpace<DISTMATDTYPE>& wspace(ReplicatedWorkSpace<DISTMATDTYPE>::instance());
    DISTMATDTYPE* work_matrix = wspace.square_matrix();

    // build a local complete matrix from a distributed matrix
    dmatrix.matgather(work_matrix, numst_);

    multiply_by_matrix(0, chromatic_number_, work_matrix, product, ldp);
}

void LocGridOrbitals::multiply_by_matrix(const int first_color,
    const int ncolors, const DISTMATDTYPE* const matrix, ORBDTYPE* product,
    const int ldp) const
{
    prod_matrix_tm_.start();

    assert(ncolors > 0);
    assert((first_color + ncolors) <= chromatic_number_);
    assert(subdivx_ > 0);

    memset(product, 0, ldp * ncolors * sizeof(ORBDTYPE));

#if 0
    (*MPIdata::sout)<<" multiply_by_matrix, first_color="<<first_color<<endl;
#endif

    DISTMATDTYPE* matrix_local = new DISTMATDTYPE[chromatic_number_ * ncolors];

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        // extract block corresponding to local indexes
        matrixToLocalMatrix(iloc, matrix, matrix_local, first_color, ncolors);

        // Compute product for subdomain iloc
        MPgemmNN(loc_numpt_, ncolors, chromatic_number_, 1., getPsi(0, iloc),
            lda_, matrix_local, chromatic_number_, 0.,
            product + iloc * loc_numpt_, ldp);
    }

#ifdef DEBUG
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {

        for (int i = 0; i < ncolors; i++)
            if (overlapping_gids_[iloc][first_color + i] == -1)
                assert(fabs(*(product + i * ldp + iloc * loc_numpt_
                              + loc_numpt_ / 2))
                       < 1.e-15);
    }
#endif

    delete[] matrix_local;

    prod_matrix_tm_.stop();
}

void LocGridOrbitals::multiplyByMatrix(const int first_color, const int ncolors,
    const SquareLocalMatrices<MATDTYPE>& matrix, ORBDTYPE* product,
    const int ldp) const
{
    prod_matrix_tm_.start();

    assert(ncolors > 0);
    assert((first_color + ncolors) <= chromatic_number_);
    assert(subdivx_ > 0);

#if 0
    (*MPIdata::sout)<<" multiplyByMatrix, first_color="<<first_color<<endl;
#endif

    // loop over subdomains
    if (chromatic_number_ > 0)
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const MATDTYPE* const mat = matrix.getSubMatrix(iloc);

            // Compute product for subdomain iloc
            MPgemmNN(loc_numpt_, chromatic_number_, ncolors, 1.,
                getPsi(0, iloc), lda_, mat + first_color, chromatic_number_, 0.,
                product + iloc * loc_numpt_, ldp);
        }

    prod_matrix_tm_.stop();
}

// Here the result is stored in one of the matrices used in the multiplication,
// so a temporary arry is necessary
void LocGridOrbitals::multiplyByMatrix(
    const SquareLocalMatrices<MATDTYPE>& matrix)
{
    prod_matrix_tm_.start();

    if (chromatic_number_ > 0)
    {
        ORBDTYPE* product = new ORBDTYPE[loc_numpt_ * chromatic_number_];
        memset(product, 0, loc_numpt_ * chromatic_number_ * sizeof(ORBDTYPE));
        const size_t slnumpt = loc_numpt_ * sizeof(ORBDTYPE);

        // loop over subdomains
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const MATDTYPE* const mat = matrix.getSubMatrix(iloc);
            ORBDTYPE* phi             = getPsi(0, iloc);

            // Compute product for subdomain iloc
            MPgemmNN(loc_numpt_, chromatic_number_, chromatic_number_, 1., phi,
                lda_, mat, chromatic_number_, 0., product, loc_numpt_);

            for (int color = 0; color < chromatic_number_; color++)
                memcpy(
                    phi + color * lda_, product + color * loc_numpt_, slnumpt);
        }

        delete[] product;
    }

    prod_matrix_tm_.stop();
}

void LocGridOrbitals::multiplyByMatrix(const int first_color, const int ncolors,
    const SquareLocalMatrices<MATDTYPE>& matrix, LocGridOrbitals& product) const
{
    multiplyByMatrix(
        first_color, ncolors, matrix, product.psi(0), product.lda_);
}

void LocGridOrbitals::multiply_by_matrix(const int first_color,
    const int ncolors, const DISTMATDTYPE* const matrix,
    LocGridOrbitals& product) const
{
    multiply_by_matrix(
        first_color, ncolors, matrix, product.psi(0), product.lda_);
}

void LocGridOrbitals::multiply_by_matrix(
    const DISTMATDTYPE* const matrix, LocGridOrbitals& product) const
{
    multiply_by_matrix(
        0, chromatic_number_, matrix, product.psi(0), product.lda_);
}

void LocGridOrbitals::multiply_by_matrix(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& matrix)
{
    prod_matrix_tm_.start();

#if 0
    (*MPIdata::sout)<<"self multiply_by_matrix"<<endl;
#endif

    ORBDTYPE* product = new ORBDTYPE[loc_numpt_ * chromatic_number_];
    memset(product, 0, loc_numpt_ * chromatic_number_ * sizeof(ORBDTYPE));

    ReplicatedWorkSpace<DISTMATDTYPE>& wspace(ReplicatedWorkSpace<DISTMATDTYPE>::instance());
    DISTMATDTYPE* work_matrix = wspace.square_matrix();

    matrix.matgather(work_matrix, numst_);

    DISTMATDTYPE* matrix_local
        = new DISTMATDTYPE[chromatic_number_ * chromatic_number_];

    const size_t slnumpt = loc_numpt_ * sizeof(ORBDTYPE);

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        ORBDTYPE* phi = getPsi(0, iloc);

        matrixToLocalMatrix(iloc, work_matrix, matrix_local);

        // Compute loc_numpt_ rows (for subdomain iloc)
        MPgemmNN(loc_numpt_, chromatic_number_, chromatic_number_, 1., phi,
            lda_, matrix_local, chromatic_number_, 0., product, loc_numpt_);

        for (int color = 0; color < chromatic_number_; color++)
            memcpy(phi + color * lda_, product + color * loc_numpt_, slnumpt);
    }

    delete[] matrix_local;
    delete[] product;

    prod_matrix_tm_.stop();
}

int LocGridOrbitals::read_hdf5(HDFrestart& h5f_file)
{
    assert(proj_matrices_ != 0);

    Control& ct = *(Control::instance());

    hid_t file_id = h5f_file.file_id();
    string name   = "Function";
    int ierr      = read_func_hdf5(h5f_file, name);
    if (ierr < 0)
    {
        (*MPIdata::serr) << "LocGridOrbitals::read_hdf5(): error in reading "
                         << name << ", size=" << name.size() << endl;
        return ierr;
    }
    else if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout) << "LocGridOrbitals::read_hdf5(): Read " << ierr
                         << " functions in restart file" << endl;
    }

    // Read DM
    if (!ct.fullyOccupied())
    {
        ierr = proj_matrices_->read_dm_hdf5(file_id);
        if (ierr < 0)
        {
            (*MPIdata::serr)
                << "LocGridOrbitals::read_hdf5(): error in reading DM" << endl;
            return ierr;
        }
    }

    return ierr;
}

int LocGridOrbitals::write_hdf5(HDFrestart& h5f_file, string name)
{
    assert(proj_matrices_ != 0);
    Control& ct = *(Control::instance());

    if (!ct.fullyOccupied())
    {
        MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
        mmpi.barrier();

        int ierr = proj_matrices_->writeDM_hdf5(h5f_file);
        if (ierr < 0) return ierr;
    }

    int ierr = write_func_hdf5(h5f_file, name);

    return ierr;
}

int LocGridOrbitals::write_func_hdf5(HDFrestart& h5f_file, string name)
{
    Control& ct   = *(Control::instance());
    hid_t file_id = h5f_file.file_id();
    bool iwrite   = h5f_file.active();

    const bool global = ct.globalColoring();
    ColoredRegions colored_regions(*pack_, *lrs_, global);

    // Create the dataspace for the dataset.

    hid_t filespace = -1;
    hid_t memspace  = -1;
    if (iwrite)
    {
        // filespace identifier
        filespace = h5f_file.createFilespace();

        // memory dataspace identifier
        memspace = h5f_file.createMemspace();
    }

    hid_t plist_id = h5f_file.createPlist();

    const short precision = ct.out_restart_info > 3 ? 2 : 1;

    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << "Write LocGridOrbitals " << name
                         << " with precision " << precision << endl;
    // loop over global (storage) functions
    for (int color = 0; color < chromatic_number_; color++)
    {
        string datasetname(getDatasetName(name, color));
        if (onpe0 && ct.verbose > 2)
            (*MPIdata::sout) << "Write " << datasetname << endl;

        // Create chunked dataset.
        hid_t dset_id = -1;

        if (iwrite)
        {
            assert(file_id > -1);

            hid_t dtype_id = ct.outHdfDataType();
            dset_id        = H5Dcreate2(file_id, datasetname.c_str(), dtype_id,
                filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
            if (dset_id < 0)
            {
                (*MPIdata::serr) << "LocGridOrbitals::write_func_hdf5(), "
                                    "H5Dcreate2 failed!!!"
                                 << endl;
                return -1;
            }
        }

        // WARNING:
        // Attributes (both the attribute information and the data it holds)
        // are considered to be metadata on an object.
        // The HDF library has a requirement that all metadata updates be done
        // collectively so all processes see the same stream of metadata
        // updates.

        // Write list of centers and radii
        // const int nst=pack_->nb_orb(color);

        vector<double> centers_and_radii;
        int nrec = 0;
        if (h5f_file.useHdf5p())
        {
            nrec = colored_regions.getAllCentersAndRadii4color(
                color, centers_and_radii);
        }
        else
        {
            nrec = colored_regions.getLocCentersAndRadii4color(
                color, centers_and_radii);
            assert(centers_and_radii.size() == 4 * nrec);
        }

        // for(int i=0;i<(int)centers_and_radii.size();i++)
        //    (*MPIdata::sout)<<"centers_and_radii["<<i<<"]="<<centers_and_radii[i]<<endl;
        // if( nrec != nst ){
        //    (*MPIdata::serr)<<"Wrong number of loc_centers in pack_!!!"<<endl;
        //}

        vector<int> gids;
        int ngids = 0;
        if (h5f_file.useHdf5p())
        {
            ngids = colored_regions.getAllGids4color(color, gids);
        }
        else
        {
            ngids = colored_regions.getLocGids4color(color, gids);
        }

        if (iwrite)
        {
            writeListCentersAndRadii(dset_id, nrec, centers_and_radii);

            writeGids(dset_id, gids);

            // Write the attribute "Lattice parameters" at "Cell origin"
            string attname("Lattice parameters");

            // Create the data space for the attribute "Lattice parameters".
            vector<double> attr_data(3);
            attr_data[0] = grid_.ll(0);
            attr_data[1] = grid_.ll(1);
            attr_data[2] = grid_.ll(2);

            mgmol_tools::addAttribute2Dataset(
                dset_id, attname.c_str(), attr_data);

            attr_data[0] = grid_.origin(0);
            attr_data[1] = grid_.origin(1);
            attr_data[2] = grid_.origin(2);

            string attname2("Cell origin");
            mgmol_tools::addAttribute2Dataset(
                dset_id, attname2.c_str(), attr_data);
        } // iwrite

        int ierr = h5f_file.writeData(
            psi(color), filespace, memspace, dset_id, precision);
        if (ierr < 0) return ierr;

        // Close/release resources.
        if (iwrite)
        {
            herr_t status = H5Dclose(dset_id);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "LocGridOrbitals::write_func_hdf5:H5Dclose failed!!!"
                    << endl;
                return -1;
            }
        }

    } // loop over color

    h5f_file.releasePlist(plist_id);

    if (iwrite)
    {
        // close filespace and memspace
        herr_t status = H5Sclose(filespace);
        if (status < 0)
        {
            (*MPIdata::serr) << "H5Sclose filespace failed!!!" << endl;
        }
        status = H5Sclose(memspace);
        if (status < 0)
        {
            (*MPIdata::serr) << "H5Sclose memspace failed!!!" << endl;
        }
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.barrier();

    return 0;
}

int LocGridOrbitals::read_func_hdf5(HDFrestart& h5f_file, string name)
{
    assert(chromatic_number_ >= 0);
    assert(name.size() > 0);
    assert(pack_);

    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    const bool global = ct.globalColoring();
    ColoredRegions colored_regions(*pack_, *lrs_, global);

    hsize_t block[3]  = { grid_.dim(0), grid_.dim(1), grid_.dim(2) };
    hsize_t offset[3] = { 0, 0, 0 };
    if (h5f_file.gatherDataX())
    {
        block[0]  = grid_.gdim(0);
        offset[1] = grid_.istart(1);
        offset[2] = grid_.istart(2);
    }

    // Each process defines dataset in memory and writes it to the hyperslab
    // in the file.

    // memory dataspace identifier
    hid_t memspace;
    if (h5f_file.active()) memspace = h5f_file.createMemspace();

    ORBDTYPE* buffer = new ORBDTYPE[block[0] * block[1] * block[2]];

    if (onpe0 && ct.verbose > 2)
    {
        if (h5f_file.gatherDataX())
        {
            (*MPIdata::sout) << "LocGridOrbitals::read_func_hdf5(): Read wave "
                                "functions from "
                             << grid_.mype_env().n_mpi_task(1)
                                    * grid_.mype_env().n_mpi_task(2)
                             << " PEs" << endl;
        }
        else
        {
            (*MPIdata::sout)
                << "LocGridOrbitals::read_func_hdf5(): Read wave functions "
                << name << " from all tasks..." << endl;
        }
    }

    vector<set<int>> filled; // set of functions already filled by data
    filled.resize(subdivx_);
    int dims[2] = { 0, 0 };

    // get centers corresponding to dataset (stored function) from input file
    multimap<string, Vector3D> centers_in_dataset;
    int ncenters = h5f_file.getLRCenters(centers_in_dataset, numst_, name);
    if (ncenters < 0) return ncenters;

    SubCell sub_cell(grid_, subdivx_, 0);
    const short precision = ct.restart_info > 3 ? 2 : 1;

    // read one color/dataset at a time
    multimap<string, Vector3D>::iterator itcenter = centers_in_dataset.begin();
    while (itcenter != centers_in_dataset.end())
    {
        vector<double> attr_data;
        short attribute_length = 0;

        const string key(itcenter->first);

        // checkif dataset exists...
        int err_id = h5f_file.dset_exists(key);
        if (h5f_file.gatherDataX()) mmpi.bcast(&err_id, 1);
        if (err_id == 0) break; // dataset does not exists

        if (onpe0 && ct.verbose > 2)
            (*MPIdata::sout) << "Read Dataset " << key << " with precision "
                             << precision << endl;

        // Open dataset.
        hid_t dset_id = h5f_file.open_dset(key);
        if (dset_id < 0)
        {
            (*MPIdata::serr)
                << "LocGridOrbitals::read_func_hdf5() --- cannot open " << key
                << endl;
            return dset_id;
        }

        herr_t status = h5f_file.readData(buffer, memspace, dset_id, precision);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "LocGridOrbitals::read_func_hdf5() --- H5Dread failed!!!"
                << endl;
            return -1;
        }

        if (h5f_file.active())
        {
            int natt = readListCentersAndRadii(dset_id, attr_data);
            assert(natt == centers_in_dataset.count(key));

            if (natt < 0) return -1;
            dims[0]          = natt;
            dims[1]          = dims[0] > 0 ? attr_data.size() / dims[0] : 0;
            attribute_length = dims[1];
        }

        status = h5f_file.close_dset(dset_id);
        if (status < 0)
        {
            return status;
        }

        int intdims[2] = { dims[0], dims[1] };

#ifdef USE_MPI
        if (h5f_file.gatherDataX())
        {
            // Bcast size of data and data
            mmpi.bcast(intdims, 2);
            const int dim = intdims[0] * intdims[1];
            assert(dim > 0);
            assert(intdims[1] == 4);
            if (!mmpi.instancePE0()) attr_data.resize(dim);
            mmpi.bcast(&attr_data[0], dim);
        }
#endif

        // loop over centers just read
        for (int i = 0; i < intdims[0]; i++)
        {
            assert(attribute_length > 0);
            assert(attr_data.size() > attribute_length * i);
            const Vector3D center(attr_data[attribute_length * i],
                attr_data[attribute_length * i + 1],
                attr_data[attribute_length * i + 2]);
            const double read_radius = attr_data[attribute_length * i + 3];

            // get possible colors for function centered at center
            set<int> possible_colors;
            colored_regions.getPossibleColors(center, possible_colors);
            // if( possible_colors.size()!=1 ){
            //    (*MPIdata::serr)<<"possible_colors.size()="<<possible_colors.size()<<endl;
            //    (*MPIdata::serr)<<"center: "<<center<<endl;
            //}
            if (possible_colors.size()
                > 0) // read center could not be local anymore...
                for (short iloc = 0; iloc < subdivx_; iloc++)
                {
                    // is this suddomain close enough to center to get data?
                    if (sub_cell.spherePossibleOverlap(
                            center, read_radius, iloc, ct.bcPoisson))
                    {
                        set<int> result;
                        // Difference of two sorted ranges
                        set_difference(possible_colors.begin(),
                            possible_colors.end(), filled[iloc].begin(),
                            filled[iloc].end(),
                            inserter(result, result.begin()));
                        if (result.size() > 0)
                        {
                            const int mycolor = (*result.begin());
                            assert(mycolor < chromatic_number_);

                            const int gid = overlapping_gids_[iloc][mycolor];
                            if (gid != -1) // gid may be -1 if no actual mesh
                                           // point is inside sphere!
                            {
                                if (masks4orbitals_->center(gid) == center)
                                {
                                    const int shift = iloc * loc_numpt_;
                                    block_vector_.assignLocal(
                                        mycolor, iloc, buffer + shift);
                                    filled[iloc].insert(mycolor);
                                    //(*MPIdata::sout)<<"Center "<<center<<",
                                    //radius "<<read_radius<<", iloc="<<iloc
                                }
                            } // gid!=-1
                            //(*MPIdata::sout)<<" Put data into mycolor
                            //"<<mycolor<<endl;
                        }
                        else
                        {
                            (*MPIdata::serr)
                                << "result.size()=" << result.size()
                                << ", possible_color="
                                << *possible_colors.begin() << "!!!" << endl;
                        }
                    } // overlap
                }

            assert(itcenter != centers_in_dataset.upper_bound(key));

            itcenter++; // increment multimap
        }
        assert(itcenter == centers_in_dataset.upper_bound(key));

    } // end loop over centers_in_dataset

    delete[] buffer;
    resetIterativeIndex();

    if (h5f_file.active())
    {
        herr_t status = H5Sclose(memspace);
        if (status < 0)
        {
            (*MPIdata::serr) << "H5Sclose failed!!!" << endl;
            return -1;
        }
    }

    return centers_in_dataset.size();
}

// initialize matrix chromatic_number_ by ncolor (for columns first_color to
// first_color+ncolor)
void LocGridOrbitals::matrixToLocalMatrix(const short iloc,
    const DISTMATDTYPE* const matrix, DISTMATDTYPE* const lmatrix) const
{
    matrixToLocalMatrix(iloc, matrix, lmatrix, 0, chromatic_number_);
}

void LocGridOrbitals::matrixToLocalMatrix(const short iloc,
    const DISTMATDTYPE* const matrix, DISTMATDTYPE* const lmatrix,
    const int first_color, const int ncolor) const
{
    assert(ncolor <= chromatic_number_);
    memset(lmatrix, 0, chromatic_number_ * ncolor * sizeof(DISTMATDTYPE));

    for (int jcolor = 0; jcolor < ncolor; jcolor++)
    {
        const int gidj = overlapping_gids_[iloc][first_color + jcolor];
        if (gidj != -1)
        {
            const int njst = gidj * numst_;
            const int njc  = jcolor * chromatic_number_;
            for (int icolor = 0; icolor < chromatic_number_; icolor++)
            {
                const int gidi = overlapping_gids_[iloc][icolor];
                if (gidi != -1)
                {
                    lmatrix[njc + icolor] = matrix[njst + gidi];
                }
            }
        }
    }
}

// compute the matrix <psi1|B|psi2>
// output: matB
void LocGridOrbitals::computeMatB(
    const LocGridOrbitals& orbitals, const pb::Lap<ORBDTYPE>& LapOper)
{
    if (numst_ == 0) return;

    assert(proj_matrices_ != 0);

    matB_tm_.start();
#if DEBUG
    if (onpe0) (*MPIdata::sout) << "LocGridOrbitals::computeMatB()" << endl;
#endif

    const short bcolor = 32;

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, chromatic_number_);

    ORBDTYPE* work = new ORBDTYPE[lda_ * bcolor];
    memset(work, 0, lda_ * bcolor * sizeof(ORBDTYPE));

    const ORBDTYPE* const orbitals_psi
        = (chromatic_number_ > 0) ? orbitals.block_vector_.vect(0) : 0;

    setDataWithGhosts();
    trade_boundaries();

    for (int icolor = 0; icolor < chromatic_number_; icolor += bcolor)
    {
        int nf = bcolor;
        if ((icolor + nf) > chromatic_number_) nf = chromatic_number_ - icolor;

        // Compute nf columns of B|psi> and store it into work
        for (int i = 0; i < nf; i++)
        {
            LapOper.rhs(getFuncWithGhosts(icolor + i), work + i * lda_);
        }

        for (short iloc = 0; iloc < subdivx_; iloc++)
        {

            MATDTYPE* ssiloc = ss.getSubMatrix(iloc);

            // calculate nf columns of ssiloc
            MPgemmTN(chromatic_number_, nf, loc_numpt_, 1.,
                orbitals_psi + iloc * loc_numpt_, lda_,
                work + iloc * loc_numpt_, lda_, 0.,
                ssiloc + icolor * chromatic_number_, chromatic_number_);
        }
    }

    delete[] work;

    const double vel = grid_.vel();
    ss.scal(vel);
    proj_matrices_->initializeMatB(ss);

    matB_tm_.stop();
}

// compute <Phi|B|Phi> and its inverse
void LocGridOrbitals::computeBAndInvB(const pb::Lap<ORBDTYPE>& LapOper)
{
    assert(proj_matrices_ != 0);

    Control& ct = *(Control::instance());
    if (!ct.Mehrstellen()) return;

    invBmat_tm_.start();

    computeMatB(*this, LapOper);
    proj_matrices_->computeInvB();

    invBmat_tm_.stop();
}

void LocGridOrbitals::getLocalOverlap(SquareLocalMatrices<MATDTYPE>& ss)
{
    assert(chromatic_number_ >= 0);
    assert(loc_numpt_ > 0);
    assert(grid_.vel() > 1.e-8);
    assert(subdivx_ > 0);

    if (chromatic_number_ != 0)
    {
#ifdef USE_MP
        getLocalOverlap(*this, ss);
#else
        const ORBDTYPE* const psi = block_vector_.vect(0);

        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            ss.syrk(iloc, loc_numpt_, psi + iloc * loc_numpt_, lda_);
        }

        // We may need the full matrix
        ss.fillUpperWithLower();

        ss.scal(grid_.vel());
#endif
    }
}

void LocGridOrbitals::getLocalOverlap(
    const LocGridOrbitals& orbitals, SquareLocalMatrices<MATDTYPE>& ss)
{
    assert(chromatic_number_ >= 0);

    if (chromatic_number_ != 0)
    {
        computeLocalProduct(
            orbitals.block_vector_.vect(0), orbitals.lda_, ss, false);
    }
}

void LocGridOrbitals::computeLocalProduct(const LocGridOrbitals& orbitals,
    LocalMatrices<MATDTYPE>& ss, const bool transpose)
{
    // assert( orbitals.chromatic_number_>=0 );
    assert(orbitals.lda_ > 1);

    if (chromatic_number_ != 0)
        computeLocalProduct(orbitals.psi(0), orbitals.lda_, ss, transpose);
}

void LocGridOrbitals::computeLocalProduct(const ORBDTYPE* const array,
    const int ld, LocalMatrices<MATDTYPE>& ss, const bool transpose)
{
    assert(loc_numpt_ > 0);
    assert(loc_numpt_ <= ld);
    assert(array != 0);
    assert(chromatic_number_ != 0);
    assert(grid_.vel() > 0.);
    assert(subdivx_ > 0);

    const ORBDTYPE* const a = transpose ? array : block_vector_.vect(0);
    const ORBDTYPE* const b = transpose ? block_vector_.vect(0) : array;

    const int lda = transpose ? ld : lda_;
    const int ldb = transpose ? lda_ : ld;

#ifdef USE_MP
    // use temporary float data for matrix ss
    LocalMatrices<ORBDTYPE> ssf(ss.subdiv(), ss.m(), ss.n());
#else
    LocalMatrices<ORBDTYPE>& ssf(ss);
#endif
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        ssf.gemm(iloc, loc_numpt_, a + iloc * loc_numpt_, lda,
            b + iloc * loc_numpt_, ldb);
    }
#ifdef USE_MP
    ss.copy(ssf);
#endif

    ss.scal(grid_.vel());
}

void LocGridOrbitals::computeDiagonalElementsDotProduct(
    const LocGridOrbitals& orbitals, vector<DISTMATDTYPE>& ss)
{
    assert(numst_ > 0);
    assert(grid_.vel() > 0.);

    memset(&ss[0], 0, numst_ * sizeof(DISTMATDTYPE));

    for (int icolor = 0; icolor < chromatic_number_; icolor++)
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const int gid = overlapping_gids_[iloc][icolor];
            if (gid > -1)
            {
                double alpha = MPdot(loc_numpt_, orbitals.getPsi(icolor, iloc),
                    getPsi(icolor, iloc));

                ss[gid] += (DISTMATDTYPE)(alpha * grid_.vel());
            }
        }
    vector<DISTMATDTYPE> tmp(ss);
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&tmp[0], &ss[0], numst_, MPI_SUM);
}

void LocGridOrbitals::computeDiagonalElementsDotProductLocal(
    const LocGridOrbitals& orbitals, vector<DISTMATDTYPE>& ss)
{
    assert(grid_.vel() > 0.);

    /* get locally centered functions */
    std::vector<int> locfcns;
    if (local_cluster_ != 0)
    {
        locfcns = local_cluster_->getClusterIndices();
    }
    else
    {
        locfcns = local_gids_;
    }

    /* define data distribution variables */
    const int siz = (int)locfcns.size();
    VariableSizeMatrix<sparserow> diag("diagDot", chromatic_number_ * subdivx_);
    /* initialize sparse diagonal matrix with locally centered indixes */
    diag.setupSparseRows(locfcns);

    /* assemble diagonal matrix */
    //    int ione=1;
    for (int icolor = 0; icolor < chromatic_number_; icolor++)
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const int ifunc = overlapping_gids_[iloc][icolor];
            if (ifunc > -1)
            {
                double alpha = MPdot(loc_numpt_, orbitals.getPsi(icolor, iloc),
                    getPsi(icolor, iloc));

                double val = alpha * grid_.vel();
                diag.insertMatrixElement(ifunc, ifunc, val, ADD, true);
            }
        }
#ifdef USE_MPI
    /* do data distribution to update local sums */
    (*distributor_diagdotprod_).augmentLocalData(diag, false);
    /* collect data */
    ss.clear();
    for (int row = 0; row < siz; row++)
    {
        ss.push_back((DISTMATDTYPE)diag.getRowEntry(row, 0));
    }
#endif
}

void LocGridOrbitals::computeGram(
    dist_matrix::DistMatrix<DISTMATDTYPE>& gram_mat)
{
    assert(proj_matrices_ != 0);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, chromatic_number_);

    getLocalOverlap(ss);

    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    gram_mat = projmatrices->getDistMatrixFromLocalMatrices(ss);
}

void LocGridOrbitals::computeGram(const LocGridOrbitals& orbitals,
    dist_matrix::DistMatrix<DISTMATDTYPE>& gram_mat)
{
    assert(proj_matrices_ != 0);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, chromatic_number_);

    getLocalOverlap(orbitals, ss);

    // make a DistMatrix out of ss
    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    gram_mat = projmatrices->getDistMatrixFromLocalMatrices(ss);
}

// compute the lower-triangular part of the overlap matrix
void LocGridOrbitals::computeGram(const int verbosity)
{
    assert(proj_matrices_ != 0);

    // if( chromatic_number_==0 )return;

    overlap_tm_.start();

#ifdef PRINT_OPERATIONS
    if (onpe0) (*MPIdata::sout) << "LocGridOrbitals::computeGram()" << endl;
#endif

    assert(subdivx_ > 0);
    assert(subdivx_ < 1000);
    assert(chromatic_number_ >= 0);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, chromatic_number_);

    getLocalOverlap(ss);

    proj_matrices_->initializeGramMatrix(ss, getIterativeIndex());

    if (verbosity > 1) proj_matrices_->printS((*MPIdata::sout));

    overlap_tm_.stop();
}

void LocGridOrbitals::computeGramAndInvS(const int verbosity)
{
    assert(proj_matrices_ != 0);

    computeGram(verbosity);

    /* Compute inverse of Gram matrix */
    proj_matrices_->computeInvS();
}

void LocGridOrbitals::checkCond(const double tol, const bool flag_stop)
{
    assert(proj_matrices_ != 0);

    proj_matrices_->checkCond(tol, flag_stop);
}

double LocGridOrbitals::dotProductWithDM(const LocGridOrbitals& orbitals)
{
    assert(proj_matrices_ != 0);
    assert(chromatic_number_ == orbitals.chromatic_number_);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, chromatic_number_);

    computeLocalProduct(orbitals, ss);

    return proj_matrices_->dotProductWithDM(ss);
}

double LocGridOrbitals::dotProductWithInvS(const LocGridOrbitals& orbitals)
{
    assert(proj_matrices_ != 0);
    assert(chromatic_number_ == orbitals.chromatic_number_);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, chromatic_number_);

    computeLocalProduct(orbitals, ss);

    return proj_matrices_->dotProductWithInvS(ss);
}

double LocGridOrbitals::dotProductDiagonal(const LocGridOrbitals& orbitals)
{
    assert(proj_matrices_ != 0);
    //(*MPIdata::sout)<<"call LocGridOrbitals::dotProductDiagonal()..."<<endl;

    vector<DISTMATDTYPE> ss;
    Control& ct = *(Control::instance());
    if (ct.short_sighted)
    {
        computeDiagonalElementsDotProductLocal(orbitals, ss);
    }
    else
    {
        ss.resize(numst_);
        computeDiagonalElementsDotProduct(orbitals, ss);
    }
    return proj_matrices_->getTraceDiagProductWithInvS(ss);
}

double LocGridOrbitals::dotProductSimple(const LocGridOrbitals& orbitals)
{
    assert(proj_matrices_ != 0);
    assert(chromatic_number_ == orbitals.chromatic_number_);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, chromatic_number_);

    computeLocalProduct(orbitals, ss);

    return proj_matrices_->dotProductSimple(ss);
}

double LocGridOrbitals::dotProduct(const LocGridOrbitals& orbitals)
{
    return (this->*dotProduct_)(orbitals); // call through pointer member
}

double LocGridOrbitals::dotProduct(
    const LocGridOrbitals& orbitals, const short dot_type)
{
    dot_product_tm_.start();

    assert(chromatic_number_ >= 0);
    assert(subdivx_ > 0);
    assert(subdivx_ < 1000);

    double dot = 0.;
    if (dot_type == 0)
    {
        dot = dotProductDiagonal(orbitals);
    }
    else if (dot_type == 1)
    {
        dot = dotProductWithInvS(orbitals);
    }
    else if (dot_type == 2)
    {
        dot = dotProductWithDM(orbitals);
    }
    else if (dot_type == 3)
    {
        dot = dotProductSimple(orbitals);
    }
    else
    {
        (*MPIdata::serr)
            << "LocGridOrbitals::dot_product() --- unknown dot product type"
            << endl;
        Control& ct = *(Control::instance());
        ct.global_exit(2);
    }

    dot_product_tm_.stop();

    return dot;
}

const dist_matrix::DistMatrix<DISTMATDTYPE> LocGridOrbitals::product(
    const LocGridOrbitals& orbitals, const bool transpose)
{
    assert(chromatic_number_ > 0);
    assert(subdivx_ > 0);
    assert(subdivx_ < 1000);

    return product(
        orbitals.psi(0), orbitals.chromatic_number_, orbitals.lda_, transpose);
}

const dist_matrix::DistMatrix<DISTMATDTYPE> LocGridOrbitals::product(
    const ORBDTYPE* const array, const int ncol, const int lda,
    const bool transpose)
{
    assert(lda > 1);

    dot_product_tm_.start();

    LocalMatrices<MATDTYPE> ss(subdivx_, chromatic_number_, ncol);

    if (chromatic_number_ != 0) computeLocalProduct(array, lda, ss, transpose);

    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    assert(projmatrices);
    dist_matrix::DistMatrix<DISTMATDTYPE> tmp(
        projmatrices->getDistMatrixFromLocalMatrices(ss));

    dot_product_tm_.stop();

    return tmp;
}

void LocGridOrbitals::orthonormalize(const bool overlap_uptodate)
{
    // if( chromatic_number_==0 )return;
    Control& ct = *(Control::instance());

    if (!overlap_uptodate) computeGram(0);

    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    assert(projmatrices);

    projmatrices->updateSubMatLS();

    const dist_matrix::SubMatrices<DISTMATDTYPE>& submatLS(
        projmatrices->getSubMatLS());
    // submatLS->print((*MPIdata::sout));

    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "LocGridOrbitals::orthonormalize()" << endl;

    for (short iloc = 0; iloc < subdivx_; iloc++)
    {

        // Loop over the functions
        for (int jcolor = 0; jcolor < chromatic_number_; jcolor++)
        {

            // compute non-diagonal elements
            for (int icolor = 0; icolor < jcolor; icolor++)
            {
                double beta
                    = (double)(-1.
                               * submatLS.val(chromatic_number_ - jcolor - 1,
                                     chromatic_number_ - icolor - 1, iloc));
                //(*MPIdata::sout)<<"beta="<<beta<<endl;
                block_vector_.axpy(beta, chromatic_number_ - icolor - 1,
                    chromatic_number_ - jcolor - 1, iloc);
            }

            // normalize state
            double alpha
                = (double)(1.
                           / submatLS.val(chromatic_number_ - jcolor - 1,
                                 chromatic_number_ - jcolor - 1, iloc));
            //(*MPIdata::sout)<<"alpha="<<alpha<<endl;
            block_vector_.scal(alpha, chromatic_number_ - jcolor - 1, iloc);
        }
    }

    incrementIterativeIndex();

    projmatrices->setGram2Id(getIterativeIndex());
}

void LocGridOrbitals::orthonormalizeLoewdin(
    const bool overlap_uptodate, SquareLocalMatrices<MATDTYPE>* matrixTransform)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "LocGridOrbitals::orthonormalizeLoewdin()" << endl;

    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    assert(projmatrices);

    if (!overlap_uptodate) computeGram(0);

    SquareLocalMatrices<MATDTYPE>* localP = matrixTransform;
    if (matrixTransform == 0)
        localP = new SquareLocalMatrices<MATDTYPE>(subdivx_, chromatic_number_);

    projmatrices->getLoewdinTransform(*localP);

    multiplyByMatrix(*localP);

    incrementIterativeIndex();

#if 0 // test
    computeGram(0);
    if( onpe0 && ct.verbose>2 )
        (*MPIdata::sout)<<"LocGridOrbitals::orthonormalizeLoewdin() --- Gram matrix (after):"<<endl;
    proj_matrices_->printS(*MPIdata::sout);
#endif
    projmatrices->setGram2Id(getIterativeIndex());

    if (matrixTransform == 0) delete localP;
}

double LocGridOrbitals::norm() const
{
    Control& ct = *(Control::instance());

    double norm = 0;

    for (int gid = 0; gid < ct.numst; gid++)
    {
        norm += normState(gid);
    }
    return norm;
}

double LocGridOrbitals::normState(const int gid) const
{
    assert(gid >= 0);

    // find color of state
    int color_gid = -1;
    for (short iloc = 0; iloc < subdivx_; iloc++)
        for (int color = 0; color < chromatic_number_; color++)
        {
            const int pst = overlapping_gids_[iloc][color];
            if (pst == gid)
            {
                assert(color_gid == -1 || color_gid == color);
                color_gid = color;
            }
        }
    double tmp = 0.;
    if (color_gid >= 0)
        for (short iloc = 0; iloc < subdivx_; iloc++)
            if (overlapping_gids_[iloc][color_gid] == gid)
            {
                // diagonal element
                tmp += block_vector_.dot(color_gid, color_gid, iloc);
                // cout<<"gid="<<gid<<", tmp="<<tmp<<endl;
            }

    double norm     = 0.;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&tmp, &norm, 1, MPI_SUM);

    return grid_.vel() * norm;
}

void LocGridOrbitals::orthonormalize2states(const int st1, const int st2)
{
    assert(st1 >= 0);
    assert(st2 >= 0);

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << "LocGridOrbitals::orthonormalize2states(): " << st1
                         << " and " << st2 << endl;
    const int st[2] = { st1, st2 };

    // find color of 2 states
    int color_st[2] = { -1, -1 };
    for (short ic = 0; ic < 2; ++ic)
    {
        color_st[ic] = pack_->getColor(st[ic]);
    }

    // work to do only if one of the states exists locally
    if (color_st[0] >= 0 || color_st[1] >= 0)
    {
        // apply mask of other function
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const int stc1
                = color_st[0] >= 0 ? overlapping_gids_[iloc][color_st[0]] : -1;
            const int stc2
                = color_st[1] >= 0 ? overlapping_gids_[iloc][color_st[1]] : -1;
            if (stc1 == st1 && stc2 == st2)
            {
                (masks4orbitals_->getMask(st1))
                    .apply(psi(color_st[1]), 0, iloc, false);
                (masks4orbitals_->getMask(st2))
                    .apply(psi(color_st[0]), 0, iloc, false);
            }
            else
            {
                if (stc1 == st1 && stc2 == -1)
                    block_vector_.set_zero(color_st[0], iloc);
                if (stc1 == -1 && stc2 == st2)
                    block_vector_.set_zero(color_st[1], iloc);
            }
        }
    }

    double tmp[3]    = { 0., 0., 0. };
    const double vel = grid_.vel();

    for (short iloc = 0; iloc < subdivx_; iloc++)
        for (int ic = 0; ic < 2; ic++)
        {
            const int color_ic = color_st[ic];

            if (color_ic >= 0)
                if (overlapping_gids_[iloc][color_ic] == st[ic])
                {
                    // diagonal element
                    tmp[2 * ic]
                        += vel * block_vector_.dot(color_ic, color_ic, iloc);

                    if (ic == 1)
                    {
                        const int color_jc = color_st[0];
                        if (color_jc >= 0)
                            if (overlapping_gids_[iloc][color_jc] == st[0])
                            {
                                tmp[1] += vel
                                          * block_vector_.dot(
                                                color_ic, color_jc, iloc);
                            }
                    }
                }
        }

    double overlap[3] = { 0., 0., 0. };
    MGmol_MPI& mmpi   = *(MGmol_MPI::instance());
    mmpi.allreduce(&tmp[0], &overlap[0], 3, MPI_SUM);

    // orthogonalize second state
    double alpha = -overlap[1] / overlap[0];
    if (color_st[0] >= 0 && color_st[1] >= 0)
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            if (overlapping_gids_[iloc][color_st[0]] == st1
                && overlapping_gids_[iloc][color_st[1]] == st2)
                block_vector_.axpy(alpha, color_st[0], color_st[1], iloc);
        }

    // normalize both states
    const double alpha1 = 1. / sqrt(overlap[0]);
    const double alpha2
        = 1. / sqrt(overlap[2] - overlap[1] * overlap[1] / overlap[0]);
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        if (color_st[0] >= 0)
            if (overlapping_gids_[iloc][color_st[0]] == st1)
                block_vector_.scal(alpha1, color_st[0], iloc);

        if (color_st[1] >= 0)
            if (overlapping_gids_[iloc][color_st[1]] == st2)
                block_vector_.scal(alpha2, color_st[1], iloc);
    }

#if 1 // testing orthonormality
    tmp[0] = 0.;
    tmp[1] = 0.;
    tmp[2] = 0.;
    for (short iloc = 0; iloc < subdivx_; iloc++)
        for (int ic = 0; ic < 2; ic++)
        {
            const int color_ic = color_st[ic];

            if (color_ic >= 0)
                if (overlapping_gids_[iloc][color_ic] == st1
                    || overlapping_gids_[iloc][color_ic] == st2)
                {
                    // diagonal element
                    tmp[2 * ic]
                        += vel * block_vector_.dot(color_ic, color_ic, iloc);

                    if (ic == 1)
                    {
                        const int color_jc = color_st[0];
                        if (color_jc >= 0)
                            if (overlapping_gids_[iloc][color_jc] == st1)
                            {
                                tmp[1] += vel
                                          * block_vector_.dot(
                                                color_ic, color_jc, iloc);
                            }
                    }
                }
        }

    mmpi.allreduce(&tmp[0], &overlap[0], 3, MPI_SUM);
    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << "Gram matrix = " << overlap[0] << "," << overlap[1]
                         << "," << overlap[2] << endl;
#endif
}

void LocGridOrbitals::multiplyByMatrix2states(
    const int st1, const int st2, const double* mat, LocGridOrbitals& product)
{
    assert(st1 >= 0);
    assert(st2 >= 0);

    if (chromatic_number_ == 0) return;

    // if( onpe0 && ct.verbose>2 )
    //  (*MPIdata::sout)<<"LocGridOrbitals::multiplyByMatrix2states()"<<endl;

    int color_st1 = -1;
    int color_st2 = -1;
    for (short iloc = 0; iloc < subdivx_; iloc++)
        for (int color = 0; color < chromatic_number_; color++)
        {
            const int st = overlapping_gids_[iloc][color];
            if (st == st1)
            {
                color_st1 = color;
                product.block_vector_.set_zero(color, iloc);
            }
            if (st == st2)
            {
                color_st2 = color;
                product.block_vector_.set_zero(color, iloc);
            }
        }

    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        if (color_st1 >= 0)
            if (overlapping_gids_[iloc][color_st1] == st1)
                product.block_vector_.axpy(
                    mat[0], block_vector_, color_st1, color_st1, iloc);
        if (color_st2 >= 0)
            if (overlapping_gids_[iloc][color_st2] == st2)
                product.block_vector_.axpy(
                    mat[3], block_vector_, color_st2, color_st2, iloc);
        if (color_st1 >= 0 && color_st2 >= 0)
            if (overlapping_gids_[iloc][color_st1] == st1
                && overlapping_gids_[iloc][color_st2] == st2)
            {
                product.block_vector_.axpy(
                    mat[2], block_vector_, color_st1, color_st2, iloc);
                product.block_vector_.axpy(
                    mat[1], block_vector_, color_st2, color_st1, iloc);
            }
    }
}

void LocGridOrbitals::computeDiagonalGram(
    VariableSizeMatrix<sparserow>& diagS) const
{
    const double vel = grid_.vel();

    for (int color = 0; color < chromatic_number_; color++)
    {
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const int gid = overlapping_gids_[iloc][color];
            if (gid != -1)
            {
                double val
                    = vel * (double)block_vector_.dot(color, color, iloc);
                diagS.insertMatrixElement(gid, gid, val, ADD, true);
            }
        }
    }
    // do data distribution to update local data.
    // All PE's need to know full diagonal entries of
    // overlapping functions, hence append=true.
    distributor_normalize_->augmentLocalData(diagS, true);
#ifdef DEBUG
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
        for (int i = 0; i < lsize; i++)
            (*MPIdata::sout)
                << "i=" << i << ", diagS[i]=" << diagS.getRowEntry(i, 0)
                << endl;
#endif
}

void LocGridOrbitals::computeInvNorms2(vector<vector<double>>& inv_norms2) const
{
    const int initTabSize = 4096;
    VariableSizeMatrix<sparserow> diagS("diagS", initTabSize);

    computeDiagonalGram(diagS);

    // assign return values
    inv_norms2.resize(subdivx_);
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        inv_norms2[iloc].resize(chromatic_number_);
        for (short color = 0; color < chromatic_number_; color++)
        {
            const int gid = overlapping_gids_[iloc][color];
            if (gid != -1)
            {
                inv_norms2[iloc][color] = diagS.get_value(gid, gid);
            }
        }
    }

    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        for (short color = 0; color < chromatic_number_; color++)
        {
            inv_norms2[iloc][color] = 1. / inv_norms2[iloc][color];
        }
    }
}

void LocGridOrbitals::normalize()
{
    normalize_tm_.start();

    assert(grid_.vel() > 1.e-8);
    assert(chromatic_number_ >= 0);

    // if( onpe0 && ct.verbose>2 )
    //        (*MPIdata::sout)<<"Normalize LocGridOrbitals"<<endl;

    Control& ct = *(Control::instance());

    if (ct.short_sighted)
    {
        const int initTabSize = 4096;
        VariableSizeMatrix<sparserow> diagS("diagSforN", initTabSize);

        computeDiagonalGram(diagS);

        /* Do normalization */
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            // Loop over the functions
            for (int color = 0; color < chromatic_number_; color++)
            {
                const int gid = overlapping_gids_[iloc][color];
                if (gid != -1)
                {
                    // normalize state
                    double alpha = 1. / sqrt(diagS.get_value(gid, gid));
                    //(*MPIdata::sout)<<"alpha="<<alpha<<endl;
                    block_vector_.scal(alpha, color, iloc);
                }
            }
        }
    }
    else
    {
        const double vel = grid_.vel();
        vector<double> diagS(numst_, 0.);
        for (int color = 0; color < chromatic_number_; color++)
        {
            for (short iloc = 0; iloc < subdivx_; iloc++)
            {

                const int gid = overlapping_gids_[iloc][color];
                if (gid != -1)
                {
                    diagS[gid]
                        += vel * (double)block_vector_.dot(color, color, iloc);
                }
            }
        }

        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.split_allreduce_sums_double(&diagS[0], numst_);
#ifdef DEBUG
        if (onpe0 && ct.verbose > 2)
            for (int i = 0; i < numst_; i++)
                (*MPIdata::sout)
                    << "i=" << i << ", diagS[i]=" << diagS[i] << endl;
#endif
        for (int i = 0; i < numst_; i++)
        {
            assert(diagS[i] > 1.e-15);
            diagS[i] = 1. / sqrt(diagS[i]);
        }
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {

            // Loop over the functions
            for (int color = 0; color < chromatic_number_; color++)
            {

                const int gid = overlapping_gids_[iloc][color];

                if (gid != -1)
                {
                    // normalize state
                    double alpha = diagS[gid];
                    //(*MPIdata::sout)<<"alpha="<<alpha<<endl;
                    block_vector_.scal(alpha, color, iloc);
                }
            }
        }
    }

    incrementIterativeIndex();

    normalize_tm_.stop();
}

// modify argument orbitals, by projecting out its component
// along LocGridOrbitals
void LocGridOrbitals::projectOut(LocGridOrbitals& orbitals, const double scale)
{
    projectOut(orbitals.psi(0), lda_, scale);

#if 0
    // test if projection is now 0
    dist_matrix::DistMatrix<DISTMATDTYPE> tmatrix(product(orbitals));
    if( onpe0 )
        (*MPIdata::sout)<<"LocGridOrbitals::projectOut(), Product after projection:"<<endl;
    tmatrix.print((*MPIdata::sout),0,0,5,5);
#endif

    orbitals.incrementIterativeIndex();
}

void LocGridOrbitals::projectOut(
    ORBDTYPE* const array, const int lda, const double scale)
{
    assert(lda > 1);
    assert(loc_numpt_ > 0);
    assert(chromatic_number_ >= 0);
    assert(lda_ > loc_numpt_);

    SquareLocalMatrices<MATDTYPE> pmatrix(subdivx_, chromatic_number_);

    if (chromatic_number_ != 0) computeLocalProduct(array, lda, pmatrix, false);

        //    pmatrix.scal(grid_.vel());

#ifdef DEBUG
    (*MPIdata::sout) << "LocGridOrbitals::projectOut()" << endl;
    (*MPIdata::sout) << "Product before projection" << endl;
    pmatrix.print((*MPIdata::sout));
#endif
    proj_matrices_->applyInvS(pmatrix);

    ORBDTYPE* tproduct = new ORBDTYPE[loc_numpt_ * chromatic_number_];
    memset(tproduct, 0, loc_numpt_ * chromatic_number_ * sizeof(ORBDTYPE));

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        ORBDTYPE* phi    = getPsi(0, iloc);
        ORBDTYPE* parray = array + iloc * loc_numpt_;

        MATDTYPE* localMat_iloc = pmatrix.getSubMatrix(iloc);

        // Compute loc_numpt_ rows (for subdomain iloc)
        MPgemmNN(loc_numpt_, chromatic_number_, chromatic_number_, 1., phi,
            lda_, localMat_iloc, chromatic_number_, 0., tproduct, loc_numpt_);

        double minus = -1. * scale;
        for (int j = 0; j < chromatic_number_; j++)
            MPaxpy(
                loc_numpt_, minus, tproduct + j * loc_numpt_, parray + j * lda);
    }

    delete[] tproduct;
}

void LocGridOrbitals::initRand()
{
    Control& ct = *(Control::instance());

    const unsigned dim[3] = { grid_.dim(0), grid_.dim(1), grid_.dim(2) };

    double* xrand = new double[grid_.gdim(0)];
    double* yrand = new double[grid_.gdim(1)];
    double* zrand = new double[grid_.gdim(2)];

    const int loc_length = dim[0] / subdivx_;
    assert(loc_length > 0);
    assert(loc_length <= dim[0]);

    const int xoff = grid_.istart(0);
    const int yoff = grid_.istart(1);
    const int zoff = grid_.istart(2);

    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << " Initialize " << chromatic_number_
                         << " random global functions" << endl;

    ran0();

    // set_zero();

    const int incx = dim[1] * dim[2];
    const int incy = dim[2];

    for (int istate = 0; istate < numst_; istate++)
    {
        // Generate x, y, z random number sequences
        for (int idx = 0; idx < grid_.gdim(0); idx++)
            xrand[idx] = ran0() - 0.5;
        for (int idx = 0; idx < grid_.gdim(1); idx++)
            yrand[idx] = ran0() - 0.5;
        for (int idx = 0; idx < grid_.gdim(2); idx++)
            zrand[idx] = ran0() - 0.5;

        int n = 0;

        for (int color = 0; color < chromatic_number_; color++)
        {
            for (short iloc = 0; iloc < subdivx_; iloc++)
            {
                assert(overlapping_gids_[iloc][color] < numst_);
                const int gid = overlapping_gids_[iloc][color];
                if (gid == istate)
                {

                    for (int ix = loc_length * iloc;
                         ix < loc_length * (iloc + 1); ix++)
                        for (int iy = 0; iy < dim[1]; iy++)
                            for (int iz = 0; iz < dim[2]; iz++)
                            {
                                const double alpha = xrand[xoff + ix]
                                                     * yrand[yoff + iy]
                                                     * zrand[zoff + iz];

                                psi(color)[ix * incx + iy * incy + iz]
                                    = alpha * alpha;
                                assert((ix * incx + iy * incy + iz) < lda_);
                            }
                    n++;

#ifdef DEBUG
                }
                else if (gid == -1)
                {
                    if (!(fabs(psi(color)[iloc * loc_length * incx]) < 1.e-15))
                        (*MPIdata::sout)
                            << "color=" << color << ", iloc=" << iloc
                            << ", val=" << psi(color)[iloc * loc_length]
                            << endl;
                    assert(fabs(psi(color)[iloc * loc_length * incx]) < 1.e-15);
#endif
                }
            }
        }
#ifdef DEBUG
        // Make sure that state istate is localized somewhere
        mmpi.allreduce(&n, 1, MPI_SUM);
        if (n <= 0)
        {
            (*MPIdata::serr)
                << "ERROR: state " << istate << " not allocated" << endl;
            (*MPIdata::serr)
                << "    chromatic_number_=" << chromatic_number_ << endl;
            for (short iloc = 0; iloc < subdivx_; iloc++)
                (*MPIdata::serr)
                    << "    overlapping_gids_[" << iloc << "][" << 0
                    << "]=" << overlapping_gids_[iloc][0] << endl;
            ct.global_exit(2);
        }
#endif
    }

    delete[] xrand;
    delete[] yrand;
    delete[] zrand;

    resetIterativeIndex();
}

// Compute nstates column of Psi^T*A*Psi starting at column first_color
// WARNING: values are added to sparse_matrix!!
void LocGridOrbitals::addDotWithNcol2Matrix(const int first_color,
    const int ncolors, LocGridOrbitals& Apsi,
    dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sparse_matrix) const
{
    addDot_tm_.start();

    assert(ncolors > 0);
    assert(chromatic_number_ > 0);

#ifdef DEBUG
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << "LocGridOrbitals::addDotWithNcol2Matrix for states "
                         << first_color << " to " << first_color + ncolors - 1
                         << endl;
    for (short icolor = 0; icolor < ncolors; icolor++)
    {
        block_vector_.hasnan(icolor);
    }
    for (int icolor = 0; icolor < ncolors; icolor++)
    {
        Apsi.block_vector_.hasnan(icolor);
    }
#endif
    const double vel = grid_.vel();

    const int size_work_cols      = chromatic_number_ * ncolors;
    DISTMATDTYPE* const work_cols = new DISTMATDTYPE[size_work_cols];
    memset(work_cols, 0,
        size_work_cols * sizeof(DISTMATDTYPE)); // necessary on bgl!!

    for (short iloc = 0; iloc < subdivx_; iloc++)
    {

        MPgemmTN(chromatic_number_, ncolors, loc_numpt_, 1.,
            block_vector_.vect(0) + iloc * loc_numpt_, lda_,
            Apsi.getPsi(first_color, iloc), lda_, 0., work_cols,
            chromatic_number_);

        for (short icolor = 0; icolor < ncolors; icolor++)
        {
            const int gid1 = overlapping_gids_[iloc][icolor];
            if (gid1 != -1)
            {
                for (int jcolor = first_color;
                     jcolor < first_color + chromatic_number_; jcolor++)
                {
                    int gid2 = overlapping_gids_[iloc][jcolor];
                    if (gid2 != -1)
                    {
                        sparse_matrix.push_back(gid1, gid2,
                            work_cols[icolor + jcolor * chromatic_number_]
                                * vel);
                    }
                }
            }
        }
    }

    delete[] work_cols;

    addDot_tm_.stop();
}

void LocGridOrbitals::addDotWithNcol2Matrix(LocGridOrbitals& Apsi,
    dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sparse_matrix) const
{
    addDotWithNcol2Matrix(0, chromatic_number_, Apsi, sparse_matrix);
}

void LocGridOrbitals::computeGlobalIndexes(LocalizationRegions& lrs)
{
    all_overlapping_gids_ = lrs.getOverlapGids();

    lrs.getLocalSubdomainIndices(local_gids_);

    overlapping_gids_.clear();
    overlapping_gids_.resize(subdivx_);
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        overlapping_gids_[iloc].resize(chromatic_number_, -1);
    }
    for (vector<int>::const_iterator it = all_overlapping_gids_.begin();
         it != all_overlapping_gids_.end(); it++)
    {
        const int gid = *it;
        short color   = pack_->getColor(gid);
        assert(color < chromatic_number_);
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            if (lrs.overlapSubdiv(gid, iloc))
            {
                overlapping_gids_[iloc][color] = gid;
            }
        }
    }
}

void LocGridOrbitals::printTimers(ostream& os)
{
    matB_tm_.print(os);
    invBmat_tm_.print(os);
    overlap_tm_.print(os);
    dot_product_tm_.print(os);
    mask_tm_.print(os);
    addDot_tm_.print(os);
    prod_matrix_tm_.print(os);
    get_dm_tm_.print((*MPIdata::sout));
    assign_tm_.print((*MPIdata::sout));
    normalize_tm_.print((*MPIdata::sout));
    axpy_tm_.print((*MPIdata::sout));
}

void LocGridOrbitals::initWF(const LocalizationRegions& lrs)
{
    Control& ct = *(Control::instance());

    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << " Initialize wave functions ..." << endl;
    }
    switch (ct.init_type)
    {
        case 1:
            if (onpe0 && ct.verbose > 1)
            {
                (*MPIdata::sout) << " with Gaussian functions..." << endl;
            }
            initGauss(ct.init_rc, lrs);
            break;
        case 2:
            if (onpe0 && ct.verbose > 1)
            {
                (*MPIdata::sout) << " with Fourier basis ..." << endl;
            }
            initFourier();
            break;
        default:
            if (onpe0 && ct.verbose > 2)
            {
                (*MPIdata::sout) << " with random values ..." << endl;
            }
            initRand();

            if (ct.globalColoring())
            {
                // smooth out random functions
                pb::Laph4M<ORBDTYPE> myoper(grid_);
                pb::GridFunc<ORBDTYPE> gf_work(
                    grid_, ct.bc[0], ct.bc[1], ct.bc[2]);
                pb::GridFunc<ORBDTYPE> gf_psi(
                    grid_, ct.bc[0], ct.bc[1], ct.bc[2]);

                if (onpe0 && ct.verbose > 2)
                    (*MPIdata::sout)
                        << " Apply B to initial wave functions" << endl;
                for (short icolor = 0; icolor < chromatic_number_; icolor++)
                {
                    gf_psi.assign(psi(icolor));
                    myoper.rhs(gf_psi, gf_work);
                    setPsi(gf_work, icolor);
                    // gf_work.init_vect(psi(icolor),'d');
                }
            }
    }
    resetIterativeIndex();

    // apply masks
    if (ct.init_type == 0)
        applyMask(true);
    else
        applyCorrMask(true);

    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout)
            << " Normalize or Orthonormalize initial wave functions" << endl;
    if (ct.isLocMode())
    {
        normalize();
        // ortho_norm_local();
    }
    else
    {
        orthonormalize();
        // orthonormalizeLoewdin();
    }

    setDataWithGhosts();
    trade_boundaries();

#ifdef DEBUG
    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << "LocGridOrbitals::init_wf() done" << endl;
#endif
}

template void LocGridOrbitals::setDataWithGhosts(
    pb::GridFuncVector<float>* data_wghosts);
template void LocGridOrbitals::setDataWithGhosts(
    pb::GridFuncVector<double>* data_wghosts);

template void LocGridOrbitals::setPsi(
    const pb::GridFunc<float>& gf_work, const int ist);
template void LocGridOrbitals::setPsi(
    const pb::GridFunc<double>& gf_work, const int ist);

template void LocGridOrbitals::setPsi(const pb::GridFuncVector<float>& gf_work);
template void LocGridOrbitals::setPsi(
    const pb::GridFuncVector<double>& gf_work);
