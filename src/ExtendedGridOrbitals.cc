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

#include "Control.h"
#include "DistMatrix.h"
#include "ExtendedGridOrbitals.h"
#include "GridFunc.h"
#include "HDFrestart.h"
#include "Laph2.h"
#include "Laph4M.h"
#include "LocalizationRegions.h"
#include "MPIdata.h"
#include "Mesh.h"
#include "ProjectedMatrices.h"
#include "ReplicatedWorkSpace.h"
#include "SquareLocalMatrices.h"
#include "SubMatrices.h"
#include "hdf_tools.h"
#include "lapack_c.h"

#include <cmath>
#include <fstream>
#include <utility>
using namespace std;

#define ORBITAL_OCCUPATION 2.
string getDatasetName(const string& name, const int color);

short ExtendedGridOrbitals::subdivx_ = 0;
int ExtendedGridOrbitals::lda_       = 0;
int ExtendedGridOrbitals::numpt_     = 0;
int ExtendedGridOrbitals::loc_numpt_ = 0;
short ExtendedGridOrbitals::bc_[3]   = { 1, 1, 1 };
ExtendedGridOrbitalsPtrFunc ExtendedGridOrbitals::dotProduct_
    = &ExtendedGridOrbitals::dotProductDiagonal;
int ExtendedGridOrbitals::data_wghosts_index_ = -1;
int ExtendedGridOrbitals::numst_              = -1;
std::vector<std::vector<int>> ExtendedGridOrbitals::overlapping_gids_;
std::shared_ptr<DataDistribution> ExtendedGridOrbitals::distributor_;

Timer ExtendedGridOrbitals::get_dm_tm_("ExtendedGridOrbitals::get_dm");
Timer ExtendedGridOrbitals::matB_tm_("ExtendedGridOrbitals::matB");
Timer ExtendedGridOrbitals::invBmat_tm_("ExtendedGridOrbitals::invBmat");
Timer ExtendedGridOrbitals::overlap_tm_("ExtendedGridOrbitals::overlap");
Timer ExtendedGridOrbitals::dot_product_tm_(
    "ExtendedGridOrbitals::dot_product");
Timer ExtendedGridOrbitals::addDot_tm_("ExtendedGridOrbitals::addDot");
Timer ExtendedGridOrbitals::prod_matrix_tm_(
    "ExtendedGridOrbitals::prod_matrix");
Timer ExtendedGridOrbitals::assign_tm_("ExtendedGridOrbitals::assign");
Timer ExtendedGridOrbitals::normalize_tm_("ExtendedGridOrbitals::normalize");
Timer ExtendedGridOrbitals::axpy_tm_("ExtendedGridOrbitals::axpy");

ExtendedGridOrbitals::ExtendedGridOrbitals(std::string name,
    const pb::Grid& my_grid, const short subdivx, const int numst,
    const short bc[3], ProjectedMatricesInterface* proj_matrices,
    LocalizationRegions* lrs, MasksSet* masks, MasksSet* corrmasks,
    ClusterOrbitals* local_cluster, const bool setup_flag)
    : name_(std::move(name)),
      proj_matrices_(proj_matrices),
      local_cluster_(local_cluster),
      block_vector_(my_grid, subdivx, bc),
      grid_(my_grid),
      lrs_(lrs)
{
    (void)masks;
    (void)corrmasks;

    // preconditions
    assert(subdivx > 0);
    assert(proj_matrices != 0);

    for (short i = 0; i < 3; i++)
        assert(bc[i] == 0 || bc[i] == 1);
    assert(grid_.size() > 0);

    subdivx_   = subdivx;
    numst_     = numst;
    numpt_     = grid_.size();
    lda_       = block_vector_.getld();
    loc_numpt_ = numpt_ / subdivx_;

    assert(numst_ >= 0);

    if (setup_flag) setup();
}

ExtendedGridOrbitals::~ExtendedGridOrbitals() { assert(proj_matrices_ != 0); }

ExtendedGridOrbitals::ExtendedGridOrbitals(const std::string& name,
    const ExtendedGridOrbitals& A, const bool copy_data)
    : Orbitals(A, copy_data),
      name_(name),
      proj_matrices_(A.proj_matrices_),
      local_cluster_(A.local_cluster_),
      block_vector_(A.block_vector_, copy_data),
      grid_(A.grid_),
      lrs_(A.lrs_)
{
    // if(onpe0)cout<<"call ExtendedGridOrbitals(const ExtendedGridOrbitals &A,
    // const bool copy_data)"<<endl;

    assert(A.proj_matrices_ != 0);
}

ExtendedGridOrbitals::ExtendedGridOrbitals(const std::string& name,
    const ExtendedGridOrbitals& A, ProjectedMatricesInterface* proj_matrices,
    const bool copy_data)
    : Orbitals(A, copy_data),
      name_(name),
      proj_matrices_(proj_matrices),
      local_cluster_(A.local_cluster_),
      block_vector_(A.block_vector_, copy_data),
      grid_(A.grid_),
      lrs_(A.lrs_)
{
    assert(proj_matrices != 0);

    // setup new projected_matrices object
    Control& ct = *(Control::instance());
    proj_matrices_->setup(ct.occ_width, ct.getNel(), overlapping_gids_);
}

void ExtendedGridOrbitals::copyDataFrom(const ExtendedGridOrbitals& src)
{
    assert(proj_matrices_ != 0);

    block_vector_.copyDataFrom(src.block_vector_);

    setIterativeIndex(src);
}

void ExtendedGridOrbitals::setDotProduct(const short dot_type)
{
    if (dot_type == 0)
        dotProduct_ = &ExtendedGridOrbitals::dotProductDiagonal;
    else if (dot_type == 1)
        dotProduct_ = &ExtendedGridOrbitals::dotProductWithInvS;
    else if (dot_type == 2)
        dotProduct_ = &ExtendedGridOrbitals::dotProductWithDM;
    else if (dot_type == 3)
        dotProduct_ = &ExtendedGridOrbitals::dotProductSimple;
}

void ExtendedGridOrbitals::setup()
{
    Control& ct = *(Control::instance());

    // preconditions
    assert(proj_matrices_ != 0);

    if (ct.verbose > 0)
        printWithTimeStamp(
            "ExtendedGridOrbitals::setup()...", (*MPIdata::sout));

    computeGlobalIndexes();

    bool skinny_stencil = !ct.Mehrstellen();

    block_vector_.initialize(overlapping_gids_, skinny_stencil);

    proj_matrices_->setup(ct.occ_width, ct.getNel(), overlapping_gids_);

    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };

    const double maxr = mygrid.maxDomainSize();
    distributor_.reset(new DataDistribution(
        "DistributorExtendedGridOrbitals", maxr, myPEenv, domain));

    if (ct.verbose > 0)
        printWithTimeStamp(
            "ExtendedGridOrbitals::setup() done...", (*MPIdata::sout));
}

void ExtendedGridOrbitals::reset(
    MasksSet* masks, MasksSet* corrmasks, LocalizationRegions* lrs)
{
    (void)masks;
    (void)corrmasks;
    (void)lrs;

    // free some old data
    block_vector_.clear();
    setIterativeIndex(-10);

    // reset
    setup();
}

void ExtendedGridOrbitals::assign(const ExtendedGridOrbitals& orbitals)
{
    assign_tm_.start();

    setIterativeIndex(orbitals);

    block_vector_.copyDataFrom(orbitals.block_vector_);

    assign_tm_.stop();
}

void ExtendedGridOrbitals::axpy(
    const double alpha, const ExtendedGridOrbitals& orbitals)
{
    axpy_tm_.start();

    block_vector_.axpy(alpha, orbitals.block_vector_);

    incrementIterativeIndex();

    axpy_tm_.stop();
}

void ExtendedGridOrbitals::applyMask(const bool first_time) { return; }

void ExtendedGridOrbitals::applyCorrMask(const bool first_time) { return; }

void ExtendedGridOrbitals::app_mask(
    const int color, ORBDTYPE* u, const short level) const
{
    return;
}

void ExtendedGridOrbitals::app_mask(
    const int color, pb::GridFunc<ORBDTYPE>& gu, const short level) const
{
    return;
}

void ExtendedGridOrbitals::init2zero()
{
    ORBDTYPE* ipsi = psi(0);
    memset(ipsi, 0, numst_ * numpt_ * sizeof(ORBDTYPE));
}

void ExtendedGridOrbitals::initGauss(
    const double rc, const LocalizationRegions& lrs)
{
    assert(numst_ >= 0);
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
    for (int icolor = 0; icolor < numst_; icolor++)
    {
        ORBDTYPE* ipsi = psi(icolor);
        memset(ipsi, 0, numpt_ * sizeof(ORBDTYPE));

        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const Vector3D& center(lrs.getCenter(icolor));
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
    resetIterativeIndex();
}

void ExtendedGridOrbitals::initFourier()
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

    const int cbrtncolors = (int)ceil(cbrt(numst_));

    for (int icolor = 0; icolor < numst_; icolor++)
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
            double x = start0 + iloc * dim0 * hgrid[0];
            for (int ix = iloc * dim0; ix < (iloc + 1) * dim0; ix++)
            {
                double y = start1;

                for (int iy = 0; iy < dim1; iy++)
                {
                    double z = start2;
                    for (int iz = 0; iz < dim2; iz++)
                    {
                        ipsi[ix * incx + iy * incy + iz] = (ORBDTYPE)(
                            cos(kk[0] * x) * cos(kk[1] * y) * cos(kk[2] * z));

                        z += hgrid[2];
                    }
                    y += hgrid[1];
                }
                x += hgrid[0];
            }
        }
    }
    resetIterativeIndex();
}

void ExtendedGridOrbitals::precond_smooth(ORBDTYPE* rhs, const int ld,
    const int ifirst, const int nvect, const int npower, const double alpha)
{
    assert(ld >= static_cast<int>(grid_.size()));

    pb::Laph2<ORBDTYPE> myoper(grid_);
    pb::GridFunc<ORBDTYPE> gf_w(grid_, bc_[0], bc_[1], bc_[2]);

    for (int i = 0; i < nvect; i++)
    {
        ORBDTYPE* rhsi = rhs + ld * i;

        pb::GridFunc<ORBDTYPE> gf_sd(rhsi, grid_, bc_[0], bc_[1], bc_[2]);

        myoper.smooth(gf_sd, gf_w, alpha);

        gf_w.init_vect(rhsi, 'd');
    }
}

void ExtendedGridOrbitals::multiply_by_matrix(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& dmatrix,
    ORBDTYPE* const product, const int ldp)
{
#if 0
    (*MPIdata::sout)<<"self multiply_by_matrix"<<endl;
#endif

    ReplicatedWorkSpace<DISTMATDTYPE>& wspace(
        ReplicatedWorkSpace<DISTMATDTYPE>::instance());
    DISTMATDTYPE* work_matrix = wspace.square_matrix();

    // build a local complete matrix from a distributed matrix
    dmatrix.allgather(work_matrix, numst_);

    multiply_by_matrix(work_matrix, product, ldp);
}

void ExtendedGridOrbitals::multiply_by_matrix(
    const DISTMATDTYPE* const matrix, ORBDTYPE* product, const int ldp) const
{
    prod_matrix_tm_.start();

    assert(subdivx_ > 0);

    memset(product, 0, ldp * numst_ * sizeof(ORBDTYPE));

#if 0
    (*MPIdata::sout)<<" multiply_by_matrix, first_color="<<first_color<<endl;
#endif

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        // Compute product for subdomain iloc
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(loc_numpt_, numst_,
            numst_, 1., getPsi(0, iloc), lda_, matrix, numst_, 0.,
            product + iloc * loc_numpt_, ldp);
    }

    prod_matrix_tm_.stop();
}

void ExtendedGridOrbitals::multiplyByMatrix(
    const SquareLocalMatrices<MATDTYPE>& matrix, ORBDTYPE* product,
    const int ldp) const
{
    prod_matrix_tm_.start();

#if 0
    (*MPIdata::sout)<<" multiplyByMatrix, first_color="<<first_color<<endl;
#endif

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        const MATDTYPE* const mat = matrix.getSubMatrix(iloc);

        // Compute product for subdomain iloc
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(loc_numpt_, numst_,
            numst_, 1., getPsi(0, iloc), lda_, mat, numst_, 0.,
            product + iloc * loc_numpt_, ldp);
    }

    prod_matrix_tm_.stop();
}

// Here the result is stored in one of the matrices used in the multiplication,
// so a temporary arry is necessary
void ExtendedGridOrbitals::multiplyByMatrix(
    const SquareLocalMatrices<MATDTYPE>& matrix)
{
    prod_matrix_tm_.start();

    ORBDTYPE* product = new ORBDTYPE[loc_numpt_ * numst_];
    memset(product, 0, loc_numpt_ * numst_ * sizeof(ORBDTYPE));
    const size_t slnumpt = loc_numpt_ * sizeof(ORBDTYPE);

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        const MATDTYPE* const mat = matrix.getSubMatrix(iloc);
        ORBDTYPE* phi             = getPsi(0, iloc);

        // Compute product for subdomain iloc
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(loc_numpt_, numst_,
            numst_, 1., phi, lda_, mat, numst_, 0., product, loc_numpt_);

        for (int color = 0; color < numst_; color++)
            memcpy(phi + color * lda_, product + color * loc_numpt_, slnumpt);
    }

    delete[] product;

    prod_matrix_tm_.stop();
}

void ExtendedGridOrbitals::multiplyByMatrix(
    const SquareLocalMatrices<MATDTYPE>& matrix,
    ExtendedGridOrbitals& product) const
{
    multiplyByMatrix(matrix, product.psi(0), product.lda_);
}

void ExtendedGridOrbitals::multiply_by_matrix(
    const DISTMATDTYPE* const matrix, ExtendedGridOrbitals& product) const
{
    multiply_by_matrix(matrix, product.psi(0), product.lda_);
}

void ExtendedGridOrbitals::multiply_by_matrix(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& matrix)
{
    prod_matrix_tm_.start();

#if 0
    (*MPIdata::sout)<<"self multiply_by_matrix"<<endl;
#endif

    ORBDTYPE* product = new ORBDTYPE[loc_numpt_ * numst_];
    memset(product, 0, loc_numpt_ * numst_ * sizeof(ORBDTYPE));

    ReplicatedWorkSpace<DISTMATDTYPE>& wspace(
        ReplicatedWorkSpace<DISTMATDTYPE>::instance());
    DISTMATDTYPE* work_matrix = wspace.square_matrix();

    matrix.allgather(work_matrix, numst_);

    const size_t slnumpt = loc_numpt_ * sizeof(ORBDTYPE);

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        ORBDTYPE* phi = getPsi(0, iloc);

        // Compute loc_numpt_ rows (for subdomain iloc)
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(loc_numpt_, numst_,
            numst_, 1., phi, lda_, work_matrix, numst_, 0., product,
            loc_numpt_);

        for (int color = 0; color < numst_; color++)
            memcpy(phi + color * lda_, product + color * loc_numpt_, slnumpt);
    }

    delete[] product;

    prod_matrix_tm_.stop();
}

int ExtendedGridOrbitals::read_hdf5(HDFrestart& h5f_file)
{
    assert(proj_matrices_ != 0);

    Control& ct = *(Control::instance());

    hid_t file_id = h5f_file.file_id();
    string name   = "Function";
    int ierr      = read_func_hdf5(h5f_file, name);
    if (ierr < 0)
    {
        (*MPIdata::serr)
            << "ExtendedGridOrbitals::read_hdf5(): error in reading " << name
            << ", size=" << name.size() << endl;
        return ierr;
    }
    else if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout) << "ExtendedGridOrbitals::read_hdf5(): Read " << ierr
                         << " functions in restart file" << endl;
    }

    // Read DM
    if (!ct.fullyOccupied())
    {
        ierr = proj_matrices_->read_dm_hdf5(file_id);
        if (ierr < 0)
        {
            (*MPIdata::serr)
                << "ExtendedGridOrbitals::read_hdf5(): error in reading DM"
                << endl;
            return ierr;
        }
    }

    resetIterativeIndex();

    return ierr;
}

int ExtendedGridOrbitals::write_hdf5(HDFrestart& h5f_file, string name)
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

    int ierr = write_func_hdf5(h5f_file, std::move(name));

    return ierr;
}

int ExtendedGridOrbitals::write_func_hdf5(
    HDFrestart& h5f_file, const string& name)
{
    Control& ct   = *(Control::instance());
    hid_t file_id = h5f_file.file_id();
    bool iwrite   = h5f_file.active();

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
        (*MPIdata::sout) << "Write ExtendedGridOrbitals " << name
                         << " with precision " << precision << endl;
    // loop over global (storage) functions
    for (int color = 0; color < numst_; color++)
    {
        string datasetname(getDatasetName(name, color));
        if (onpe0 && ct.verbose > 2)
            (*MPIdata::sout) << "Write " << datasetname << endl;

        // Create chunked dataset.
        hid_t dset_id = -1;

        if (iwrite)
        {
            assert(file_id > -1);

            hid_t dtype_id = outHdfDataType(ct.out_restart_info);
            dset_id        = H5Dcreate2(file_id, datasetname.c_str(), dtype_id,
                filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
            if (dset_id < 0)
            {
                (*MPIdata::serr) << "ExtendedGridOrbitals::write_func_hdf5(), "
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

        vector<int> gids;
        gids.push_back(color);

        if (iwrite)
        {
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
                (*MPIdata::serr) << "ExtendedGridOrbitals::write_func_hdf5:"
                                    "H5Dclose failed!!!"
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

// read all the data sets with names starting with "name"
int ExtendedGridOrbitals::read_func_hdf5(
    HDFrestart& h5f_file, const string& name)
{
    assert(numst_ >= 0);
    assert(name.size() > 0);

    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    hsize_t block[3] = { grid_.dim(0), grid_.dim(1), grid_.dim(2) };
    if (h5f_file.gatherDataX())
    {
        block[0] = grid_.gdim(0);
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
            (*MPIdata::sout)
                << "ExtendedGridOrbitals::read_func_hdf5(): Read wave "
                   "functions from "
                << grid_.mype_env().n_mpi_task(1)
                       * grid_.mype_env().n_mpi_task(2)
                << " PEs" << endl;
        }
        else
        {
            (*MPIdata::sout) << "ExtendedGridOrbitals::read_func_hdf5(): Read "
                                "wave functions "
                             << name << " from all tasks..." << endl;
        }
    }

    const short precision = ct.restart_info > 3 ? 2 : 1;

    for (int icolor = 0; icolor < numst_; icolor++)
    {
        const string datasetname(getDatasetName(name, icolor));

        // check if dataset exists...
        int err_id = h5f_file.dset_exists(datasetname);
        if (h5f_file.gatherDataX()) mmpi.bcast(&err_id, 1);
        if (err_id == 0) break; // dataset does not exists

        if (onpe0 && ct.verbose > 2)
            (*MPIdata::sout) << "Read Dataset " << datasetname
                             << " with precision " << precision << endl;

        // Open dataset.
        hid_t dset_id = h5f_file.open_dset(datasetname);
        if (dset_id < 0)
        {
            (*MPIdata::serr)
                << "ExtendedGridOrbitals::read_func_hdf5() --- cannot open "
                << datasetname << endl;
            return dset_id;
        }

        herr_t status = h5f_file.readData(buffer, memspace, dset_id, precision);
        if (status < 0)
        {
            (*MPIdata::serr) << "ExtendedGridOrbitals::read_func_hdf5() --- "
                                "H5Dread failed!!!"
                             << endl;
            return -1;
        }

        status = h5f_file.close_dset(dset_id);
        if (status < 0)
        {
            return status;
        }

        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const int shift = iloc * loc_numpt_;
            block_vector_.assignLocal(icolor, iloc, buffer + shift);
        }
    }

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

    return numst_;
}

// compute the matrix <psi1|B|psi2>
// output: matB
void ExtendedGridOrbitals::computeMatB(
    const ExtendedGridOrbitals& orbitals, const pb::Lap<ORBDTYPE>& LapOper)
{
    if (numst_ == 0) return;

    assert(proj_matrices_ != 0);

    matB_tm_.start();
#if DEBUG
    if (onpe0)
        (*MPIdata::sout) << "ExtendedGridOrbitals::computeMatB()" << endl;
#endif

    const short bcolor = 32;

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    ORBDTYPE* work = new ORBDTYPE[lda_ * bcolor];
    memset(work, 0, lda_ * bcolor * sizeof(ORBDTYPE));

    const ORBDTYPE* const orbitals_psi
        = (numst_ > 0) ? orbitals.block_vector_.vect(0) : nullptr;

    setDataWithGhosts();
    trade_boundaries();

    for (int icolor = 0; icolor < numst_; icolor += bcolor)
    {
        int nf = bcolor;
        if ((icolor + nf) > numst_) nf = numst_ - icolor;

        // Compute nf columns of B|psi> and store it into work
        for (int i = 0; i < nf; i++)
        {
            LapOper.rhs(getFuncWithGhosts(icolor + i), work + i * lda_);
        }

        for (short iloc = 0; iloc < subdivx_; iloc++)
        {

            MATDTYPE* ssiloc = ss.getSubMatrix(iloc);

            // calculate nf columns of ssiloc
            MPgemmTN(numst_, nf, loc_numpt_, 1.,
                orbitals_psi + iloc * loc_numpt_, lda_,
                work + iloc * loc_numpt_, lda_, 0., ssiloc + icolor * numst_,
                numst_);
        }
    }

    delete[] work;

    const double vel = grid_.vel();
    ss.scal(vel);
    proj_matrices_->initializeMatB(ss);

    matB_tm_.stop();
}

// compute <Phi|B|Phi> and its inverse
void ExtendedGridOrbitals::computeBAndInvB(const pb::Lap<ORBDTYPE>& LapOper)
{
    assert(proj_matrices_ != 0);

    Control& ct = *(Control::instance());
    if (!ct.Mehrstellen()) return;

    invBmat_tm_.start();

    computeMatB(*this, LapOper);
    proj_matrices_->computeInvB();

    invBmat_tm_.stop();
}

void ExtendedGridOrbitals::getLocalOverlap(SquareLocalMatrices<MATDTYPE>& ss)
{
    assert(numst_ >= 0);
    assert(loc_numpt_ > 0);
    assert(grid_.vel() > 1.e-8);
    assert(subdivx_ > 0);

    if (numst_ != 0)
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

void ExtendedGridOrbitals::getLocalOverlap(
    const ExtendedGridOrbitals& orbitals, SquareLocalMatrices<MATDTYPE>& ss)
{
    assert(numst_ >= 0);

    if (numst_ != 0)
    {
        computeLocalProduct(
            orbitals.block_vector_.vect(0), orbitals.lda_, ss, false);
    }
}

void ExtendedGridOrbitals::computeLocalProduct(
    const ExtendedGridOrbitals& orbitals, LocalMatrices<MATDTYPE>& ss,
    const bool transpose)
{
    // assert( orbitals.numst_>=0 );
    assert(orbitals.lda_ > 1);

    if (numst_ != 0)
        computeLocalProduct(orbitals.psi(0), orbitals.lda_, ss, transpose);
}

void ExtendedGridOrbitals::computeLocalProduct(const ORBDTYPE* const array,
    const int ld, LocalMatrices<MATDTYPE>& ss, const bool transpose)
{
    assert(loc_numpt_ > 0);
    assert(loc_numpt_ <= ld);
    assert(array != 0);
    assert(numst_ != 0);
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

void ExtendedGridOrbitals::computeDiagonalElementsDotProduct(
    const ExtendedGridOrbitals& orbitals, vector<DISTMATDTYPE>& ss) const
{
    assert(numst_ > 0);
    assert(grid_.vel() > 0.);

    for (int icolor = 0; icolor < numst_; icolor++)
    {
        ss[icolor] = 0.;
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            double alpha = MPdot(loc_numpt_, orbitals.getPsi(icolor, iloc),
                getPsi(icolor, iloc));

            ss[icolor] += (DISTMATDTYPE)(alpha * grid_.vel());
        }
    }
    vector<DISTMATDTYPE> tmp(ss);
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&tmp[0], &ss[0], numst_, MPI_SUM);
}

void ExtendedGridOrbitals::computeGram(
    dist_matrix::DistMatrix<DISTMATDTYPE>& gram_mat)
{
    assert(proj_matrices_ != 0);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    getLocalOverlap(ss);

    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    gram_mat = projmatrices->getDistMatrixFromLocalMatrices(ss);
}

void ExtendedGridOrbitals::computeGram(const ExtendedGridOrbitals& orbitals,
    dist_matrix::DistMatrix<DISTMATDTYPE>& gram_mat)
{
    assert(proj_matrices_ != 0);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    getLocalOverlap(orbitals, ss);

    // make a DistMatrix out of ss
    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    gram_mat = projmatrices->getDistMatrixFromLocalMatrices(ss);
}

// compute the lower-triangular part of the overlap matrix
void ExtendedGridOrbitals::computeGram(const int verbosity)
{
    assert(proj_matrices_ != 0);

    overlap_tm_.start();

#ifdef PRINT_OPERATIONS
    if (onpe0)
        (*MPIdata::sout) << "ExtendedGridOrbitals::computeGram()" << endl;
#endif

    assert(subdivx_ > 0);
    assert(subdivx_ < 1000);
    assert(numst_ >= 0);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    getLocalOverlap(ss);

    proj_matrices_->initializeGramMatrix(ss, getIterativeIndex());

    if (verbosity > 1) proj_matrices_->printS((*MPIdata::sout));

    overlap_tm_.stop();
}

void ExtendedGridOrbitals::computeGramAndInvS(const int verbosity)
{
    assert(proj_matrices_ != 0);

    computeGram(verbosity);

    /* Compute inverse of Gram matrix */
    proj_matrices_->computeInvS();
}

void ExtendedGridOrbitals::checkCond(const double tol, const bool flag_stop)
{
    assert(proj_matrices_ != 0);

    proj_matrices_->checkCond(tol, flag_stop);
}

double ExtendedGridOrbitals::dotProductWithDM(
    const ExtendedGridOrbitals& orbitals)
{
    assert(proj_matrices_ != 0);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    computeLocalProduct(orbitals, ss);

    return proj_matrices_->dotProductWithDM(ss);
}

double ExtendedGridOrbitals::dotProductWithInvS(
    const ExtendedGridOrbitals& orbitals)
{
    assert(proj_matrices_ != 0);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    computeLocalProduct(orbitals, ss);

    return proj_matrices_->dotProductWithInvS(ss);
}

double ExtendedGridOrbitals::dotProductDiagonal(
    const ExtendedGridOrbitals& orbitals)
{
    assert(proj_matrices_ != 0);

    vector<DISTMATDTYPE> ss(numst_);
    computeDiagonalElementsDotProduct(orbitals, ss);

    return proj_matrices_->getTraceDiagProductWithInvS(ss);
}

double ExtendedGridOrbitals::dotProductSimple(
    const ExtendedGridOrbitals& orbitals)
{
    assert(proj_matrices_ != 0);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    computeLocalProduct(orbitals, ss);

    return proj_matrices_->dotProductSimple(ss);
}

double ExtendedGridOrbitals::dotProduct(const ExtendedGridOrbitals& orbitals)
{
    return (this->*dotProduct_)(orbitals); // call through pointer member
}

double ExtendedGridOrbitals::dotProduct(
    const ExtendedGridOrbitals& orbitals, const short dot_type)
{
    dot_product_tm_.start();

    assert(numst_ >= 0);
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
        (*MPIdata::serr) << "ExtendedGridOrbitals::dot_product() --- unknown "
                            "dot product type"
                         << endl;
        Control& ct = *(Control::instance());
        ct.global_exit(2);
    }

    dot_product_tm_.stop();

    return dot;
}

const dist_matrix::DistMatrix<DISTMATDTYPE> ExtendedGridOrbitals::product(
    const ExtendedGridOrbitals& orbitals, const bool transpose)
{
    assert(numst_ > 0);
    assert(subdivx_ > 0);
    assert(subdivx_ < 1000);

    return product(orbitals.psi(0), numst_, orbitals.lda_, transpose);
}

const dist_matrix::DistMatrix<DISTMATDTYPE> ExtendedGridOrbitals::product(
    const ORBDTYPE* const array, const int ncol, const int lda,
    const bool transpose)
{
    assert(lda > 1);

    dot_product_tm_.start();

    LocalMatrices<MATDTYPE> ss(subdivx_, numst_, ncol);

    computeLocalProduct(array, lda, ss, transpose);

    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    assert(projmatrices);
    dist_matrix::DistMatrix<DISTMATDTYPE> tmp(
        projmatrices->getDistMatrixFromLocalMatrices(ss));

    dot_product_tm_.stop();

    return tmp;
}

void ExtendedGridOrbitals::orthonormalize(const bool overlap_uptodate)
{
    Control& ct = *(Control::instance());

    if (!overlap_uptodate) computeGram(0);

    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    assert(projmatrices);

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    const dist_matrix::SubMatrices<DISTMATDTYPE>& submatLS(
        projmatrices->getSubMatLS(mmpi.commSameSpin(), overlapping_gids_));
    // submatLS.print((*MPIdata::sout));

    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "ExtendedGridOrbitals::orthonormalize()" << endl;

    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        // Loop over the functions
        for (int jcolor = 0; jcolor < numst_; jcolor++)
        {
            // compute non-diagonal elements
            for (int icolor = 0; icolor < jcolor; icolor++)
            {
                double beta = (double)(-1.
                                       * submatLS.val(numst_ - jcolor - 1,
                                             numst_ - icolor - 1, iloc));
                //(*MPIdata::sout)<<"beta="<<beta<<endl;
                block_vector_.axpy(
                    beta, numst_ - icolor - 1, numst_ - jcolor - 1, iloc);
            }

            // normalize state
            double alpha = (double)(1.
                                    / submatLS.val(numst_ - jcolor - 1,
                                          numst_ - jcolor - 1, iloc));
            //(*MPIdata::sout)<<"alpha="<<alpha<<endl;
            block_vector_.scal(alpha, numst_ - jcolor - 1, iloc);
        }
    }

    incrementIterativeIndex();

    projmatrices->setGram2Id(getIterativeIndex());
}

void ExtendedGridOrbitals::orthonormalizeLoewdin(
    const bool overlap_uptodate, SquareLocalMatrices<MATDTYPE>* matrixTransform)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "ExtendedGridOrbitals::orthonormalizeLoewdin()"
                         << endl;

    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    assert(projmatrices);

    if (!overlap_uptodate) computeGram(0);

    SquareLocalMatrices<MATDTYPE>* localP = matrixTransform;
    if (matrixTransform == nullptr)
        localP = new SquareLocalMatrices<MATDTYPE>(subdivx_, numst_);

    incrementIterativeIndex();

    projmatrices->computeLoewdinTransform(*localP, getIterativeIndex());

    multiplyByMatrix(*localP);

#if 0 // test
    computeGram(0);
    if( onpe0 && ct.verbose>2 )
        (*MPIdata::sout)<<"ExtendedGridOrbitals::orthonormalizeLoewdin() --- Gram matrix (after):"<<endl;
    proj_matrices_->printS(*MPIdata::sout);
#endif
    projmatrices->setGram2Id(getIterativeIndex());

    if (matrixTransform == nullptr) delete localP;
}

double ExtendedGridOrbitals::norm() const
{
    double norm = 0;

    for (int gid = 0; gid < numst_; gid++)
    {
        norm += normState(gid);
    }
    return norm;
}

double ExtendedGridOrbitals::normState(const int gid) const
{
    assert(gid >= 0);

    double tmp = 0.;
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        // diagonal element
        tmp += block_vector_.dot(gid, gid, iloc);
        // cout<<"gid="<<gid<<", tmp="<<tmp<<endl;
    }

    double norm     = 0.;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&tmp, &norm, 1, MPI_SUM);

    return grid_.vel() * norm;
}

void ExtendedGridOrbitals::orthonormalize2states(const int st1, const int st2)
{
    assert(st1 >= 0);
    assert(st2 >= 0);

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << "ExtendedGridOrbitals::orthonormalize2states(): "
                         << st1 << " and " << st2 << endl;
    const int st[2] = { st1, st2 };

    double tmp[3]    = { 0., 0., 0. };
    const double vel = grid_.vel();

    for (short iloc = 0; iloc < subdivx_; iloc++)
        for (int ic = 0; ic < 2; ic++)
        {
            const int color_ic = st[ic];

            // diagonal element
            tmp[2 * ic] += vel * block_vector_.dot(color_ic, color_ic, iloc);

            if (ic == 1)
            {
                const int color_jc = st[0];

                tmp[1] += vel * block_vector_.dot(color_ic, color_jc, iloc);
            }
        }

    double overlap[3] = { 0., 0., 0. };
    MGmol_MPI& mmpi   = *(MGmol_MPI::instance());
    mmpi.allreduce(&tmp[0], &overlap[0], 3, MPI_SUM);

    // orthogonalize second state
    double alpha = -overlap[1] / overlap[0];
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        block_vector_.axpy(alpha, st[0], st[1], iloc);
    }

    // normalize both states
    const double alpha1 = 1. / sqrt(overlap[0]);
    const double alpha2
        = 1. / sqrt(overlap[2] - overlap[1] * overlap[1] / overlap[0]);
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        block_vector_.scal(alpha1, st[0], iloc);
        block_vector_.scal(alpha2, st[1], iloc);
    }

#if 1 // testing orthonormality
    tmp[0] = 0.;
    tmp[1] = 0.;
    tmp[2] = 0.;
    for (short iloc = 0; iloc < subdivx_; iloc++)
        for (int ic = 0; ic < 2; ic++)
        {
            const int color_ic = st[ic];

            // diagonal element
            tmp[2 * ic] += vel * block_vector_.dot(color_ic, color_ic, iloc);

            if (ic == 1)
            {
                const int color_jc = st[0];

                tmp[1] += vel * block_vector_.dot(color_ic, color_jc, iloc);
            }
        }

    mmpi.allreduce(&tmp[0], &overlap[0], 3, MPI_SUM);
    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << "Gram matrix = " << overlap[0] << "," << overlap[1]
                         << "," << overlap[2] << endl;
#endif
}

void ExtendedGridOrbitals::multiplyByMatrix2states(const int st1, const int st2,
    const double* mat, ExtendedGridOrbitals& product)
{
    assert(st1 >= 0);
    assert(st2 >= 0);

    // if( onpe0 && ct.verbose>2 )
    //  (*MPIdata::sout)<<"ExtendedGridOrbitals::multiplyByMatrix2states()"<<endl;

    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        product.block_vector_.set_zero(st1, iloc);
        product.block_vector_.set_zero(st2, iloc);
    }

    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        product.block_vector_.axpy(mat[0], block_vector_, st1, st1, iloc);
        product.block_vector_.axpy(mat[3], block_vector_, st2, st2, iloc);
        product.block_vector_.axpy(mat[2], block_vector_, st1, st2, iloc);
        product.block_vector_.axpy(mat[1], block_vector_, st2, st1, iloc);
    }
}

void ExtendedGridOrbitals::computeInvNorms2(
    vector<vector<double>>& inv_norms2) const
{
    vector<double> diagS(numst_);

    computeDiagonalElementsDotProduct(*this, diagS);

    inv_norms2.resize(subdivx_);
    for (short iloc = 0; iloc < subdivx_; iloc++)
        inv_norms2[iloc].resize(numst_);

    for (short color = 0; color < numst_; color++)
    {
        double alpha = 1. / diagS[color];

        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            inv_norms2[iloc][color] = alpha;
        }
    }
}

void ExtendedGridOrbitals::normalize()
{
    normalize_tm_.start();

    assert(grid_.vel() > 1.e-8);
    assert(numst_ >= 0);

    // if( onpe0 && ct.verbose>2 )
    //        (*MPIdata::sout)<<"Normalize ExtendedGridOrbitals"<<endl;

    //    const double vel = grid_.vel();
    vector<double> diagS(numst_);

    computeDiagonalElementsDotProduct(*this, diagS);

    for (int color = 0; color < numst_; color++)
    {
#ifdef DEBUG
        if (onpe0 && ct.verbose > 2)
            for (int i = 0; i < numst_; i++)
                (*MPIdata::sout)
                    << "i=" << i << ", diagS[i]=" << diagS[i] << endl;
#endif
        assert(diagS[color] > 1.e-15);
        diagS[color] = 1. / sqrt(diagS[color]);

        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            block_vector_.scal(diagS[color], color, iloc);
        }
    }

    incrementIterativeIndex();

    normalize_tm_.stop();
}

// modify argument orbitals, by projecting out its component
// along ExtendedGridOrbitals
void ExtendedGridOrbitals::projectOut(
    ExtendedGridOrbitals& orbitals, const double scale)
{
    projectOut(orbitals.psi(0), lda_, scale);

#if 0
    // test if projection is now 0
    dist_matrix::DistMatrix<DISTMATDTYPE> tmatrix(product(orbitals));
    if( onpe0 )
        (*MPIdata::sout)<<"ExtendedGridOrbitals::projectOut(), Product after projection:"<<endl;
    tmatrix.print((*MPIdata::sout),0,0,5,5);
#endif

    orbitals.incrementIterativeIndex();
}

void ExtendedGridOrbitals::projectOut(
    ORBDTYPE* const array, const int lda, const double scale)
{
    assert(lda > 1);
    assert(loc_numpt_ > 0);
    assert(numst_ >= 0);
    assert(lda_ > loc_numpt_);

    SquareLocalMatrices<MATDTYPE> pmatrix(subdivx_, numst_);

    if (numst_ != 0) computeLocalProduct(array, lda, pmatrix, false);

        //    pmatrix.scal(grid_.vel());

#ifdef DEBUG
    (*MPIdata::sout) << "ExtendedGridOrbitals::projectOut()" << endl;
    (*MPIdata::sout) << "Product before projection" << endl;
    pmatrix.print((*MPIdata::sout));
#endif
    proj_matrices_->applyInvS(pmatrix);

    ORBDTYPE* tproduct = new ORBDTYPE[loc_numpt_ * numst_];
    memset(tproduct, 0, loc_numpt_ * numst_ * sizeof(ORBDTYPE));

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        ORBDTYPE* phi    = getPsi(0, iloc);
        ORBDTYPE* parray = array + iloc * loc_numpt_;

        MATDTYPE* localMat_iloc = pmatrix.getSubMatrix(iloc);

        // Compute loc_numpt_ rows (for subdomain iloc)
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(loc_numpt_, numst_,
            numst_, 1., phi, lda_, localMat_iloc, numst_, 0., tproduct,
            loc_numpt_);

        double minus = -1. * scale;
        for (int j = 0; j < numst_; j++)
            MPaxpy(
                loc_numpt_, minus, tproduct + j * loc_numpt_, parray + j * lda);
    }

    delete[] tproduct;
}

void ExtendedGridOrbitals::initRand()
{
    Control& ct = *(Control::instance());

    const unsigned dim[3] = { grid_.dim(0), grid_.dim(1), grid_.dim(2) };

    double* xrand = new double[grid_.gdim(0)];
    double* yrand = new double[grid_.gdim(1)];
    double* zrand = new double[grid_.gdim(2)];

    const int loc_length = dim[0] / subdivx_;
    assert(loc_length > 0);
    assert(static_cast<unsigned int>(loc_length) <= dim[0]);

    const int xoff = grid_.istart(0);
    const int yoff = grid_.istart(1);
    const int zoff = grid_.istart(2);

    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << " Initialize " << numst_
                         << " random global functions" << endl;

    ran0();

    // set_zero();

    const int incx = dim[1] * dim[2];
    const int incy = dim[2];

    for (int istate = 0; istate < numst_; istate++)
    {
        // Generate x, y, z random number sequences
        for (unsigned int idx = 0; idx < grid_.gdim(0); idx++)
            xrand[idx] = ran0() - 0.5;
        for (unsigned int idx = 0; idx < grid_.gdim(1); idx++)
            yrand[idx] = ran0() - 0.5;
        for (unsigned int idx = 0; idx < grid_.gdim(2); idx++)
            zrand[idx] = ran0() - 0.5;

        int n = 0;

        {
            for (short iloc = 0; iloc < subdivx_; iloc++)
            {
                for (int ix = loc_length * iloc; ix < loc_length * (iloc + 1);
                     ix++)
                    for (unsigned int iy = 0; iy < dim[1]; iy++)
                        for (unsigned int iz = 0; iz < dim[2]; iz++)
                        {
                            const double alpha = xrand[xoff + ix]
                                                 * yrand[yoff + iy]
                                                 * zrand[zoff + iz];

                            psi(istate)[ix * incx + iy * incy + iz]
                                = alpha * alpha;
                            assert((ix * incx + iy * incy + iz)
                                   < static_cast<unsigned int>(lda_));
                        }
                n++;
            }
        }
    }

    delete[] xrand;
    delete[] yrand;
    delete[] zrand;

    resetIterativeIndex();
}

// Compute nstates column of Psi^T*A*Psi starting at column first_color
// WARNING: values are added to sparse_matrix!!
void ExtendedGridOrbitals::addDotWithNcol2Matrix(const int first_color,
    const int ncolors, ExtendedGridOrbitals& Apsi,
    dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sparse_matrix) const
{
    addDot_tm_.start();

    assert(ncolors > 0);
    assert(numst_ > 0);

#ifdef DEBUG
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout)
            << "ExtendedGridOrbitals::addDotWithNcol2Matrix for states "
            << first_color << " to " << first_color + ncolors - 1 << endl;
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

    const int size_work_cols      = numst_ * ncolors;
    DISTMATDTYPE* const work_cols = new DISTMATDTYPE[size_work_cols];
    memset(work_cols, 0,
        size_work_cols * sizeof(DISTMATDTYPE)); // necessary on bgl!!

    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        MPgemmTN(numst_, ncolors, loc_numpt_, 1.,
            block_vector_.vect(0) + iloc * loc_numpt_, lda_,
            Apsi.getPsi(first_color, iloc), lda_, 0., work_cols, numst_);

        for (short icolor = 0; icolor < ncolors; icolor++)
        {
            for (int jcolor = first_color; jcolor < first_color + numst_;
                 jcolor++)
            {
                sparse_matrix.push_back(
                    icolor, jcolor, work_cols[icolor + jcolor * numst_] * vel);
            }
        }
    }

    delete[] work_cols;

    addDot_tm_.stop();
}

void ExtendedGridOrbitals::addDotWithNcol2Matrix(ExtendedGridOrbitals& Apsi,
    dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sparse_matrix) const
{
    addDotWithNcol2Matrix(0, numst_, Apsi, sparse_matrix);
}

void ExtendedGridOrbitals::computeGlobalIndexes()
{
    overlapping_gids_.clear();
    overlapping_gids_.resize(subdivx_);
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        overlapping_gids_[iloc].resize(numst_, -1);
        for (int gid = 0; gid < numst_; gid++)
        {
            overlapping_gids_[iloc][gid] = gid;
        }
    }
}

void ExtendedGridOrbitals::printTimers(ostream& os)
{
    matB_tm_.print(os);
    invBmat_tm_.print(os);
    overlap_tm_.print(os);
    dot_product_tm_.print(os);
    addDot_tm_.print(os);
    prod_matrix_tm_.print(os);
    get_dm_tm_.print(os);
    assign_tm_.print(os);
    normalize_tm_.print(os);
    axpy_tm_.print(os);
}

void ExtendedGridOrbitals::initWF(const LocalizationRegions& lrs)
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
                for (short icolor = 0; icolor < numst_; icolor++)
                {
                    gf_psi.assign(psi(icolor));
                    myoper.rhs(gf_psi, gf_work);
                    setPsi(gf_work, icolor);
                    // gf_work.init_vect(psi(icolor),'d');
                }
            }
    }
    resetIterativeIndex();

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
        // orthonormalize();
        orthonormalizeLoewdin();
    }

    setDataWithGhosts();
    trade_boundaries();

#ifdef DEBUG
    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << "ExtendedGridOrbitals::init_wf() done" << endl;
#endif
}

template void ExtendedGridOrbitals::setDataWithGhosts(
    pb::GridFuncVector<float>* data_wghosts);
template void ExtendedGridOrbitals::setDataWithGhosts(
    pb::GridFuncVector<double>* data_wghosts);

template void ExtendedGridOrbitals::setPsi(
    const pb::GridFunc<float>& gf_work, const int ist);
template void ExtendedGridOrbitals::setPsi(
    const pb::GridFunc<double>& gf_work, const int ist);

template void ExtendedGridOrbitals::setPsi(
    const pb::GridFuncVector<float>& gf_work);
template void ExtendedGridOrbitals::setPsi(
    const pb::GridFuncVector<double>& gf_work);
