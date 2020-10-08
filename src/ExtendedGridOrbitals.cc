// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "global.h"

#include <mpi.h>

#include "Control.h"
#include "DistMatrix.h"
#include "ExtendedGridOrbitals.h"
#include "GridFunc.h"
#include "HDFrestart.h"
#include "Laph2.h"
#include "Laph4M.h"
#include "LocalMatrices2DistMatrix.h"
#include "LocalizationRegions.h"
#include "MPIdata.h"
#include "Mesh.h"
#include "ProjectedMatrices.h"
#include "ReplicatedWorkSpace.h"
#include "SquareLocalMatrices.h"
#include "hdf_tools.h"
#include "lapack_c.h"
#include "memory_space.h"

#include <cmath>
#include <fstream>
#include <utility>

#define ORBITAL_OCCUPATION 2.
std::string getDatasetName(const std::string& name, const int color);

short ExtendedGridOrbitals::subdivx_ = 0;
int ExtendedGridOrbitals::lda_       = 0;
int ExtendedGridOrbitals::numpt_     = 0;
int ExtendedGridOrbitals::loc_numpt_ = 0;
ExtendedGridOrbitalsPtrFunc ExtendedGridOrbitals::dotProduct_
    = &ExtendedGridOrbitals::dotProductDiagonal;
int ExtendedGridOrbitals::data_wghosts_index_ = -1;
int ExtendedGridOrbitals::numst_              = -1;
std::vector<std::vector<int>> ExtendedGridOrbitals::overlapping_gids_;

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
    std::shared_ptr<LocalizationRegions> lrs, MasksSet* masks,
    MasksSet* corrmasks, ClusterOrbitals* local_cluster, const bool setup_flag)
    : name_(std::move(name)),
      proj_matrices_(proj_matrices),
      block_vector_(my_grid, subdivx, bc),
      grid_(my_grid)
{
    (void)lrs;
    (void)masks;
    (void)corrmasks;
    (void)local_cluster;

    // preconditions
    assert(subdivx > 0);
    assert(proj_matrices != nullptr);

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

ExtendedGridOrbitals::~ExtendedGridOrbitals()
{
    assert(proj_matrices_ != nullptr);
}

ExtendedGridOrbitals::ExtendedGridOrbitals(const std::string& name,
    const ExtendedGridOrbitals& A, const bool copy_data)
    : Orbitals(A, copy_data),
      name_(name),
      proj_matrices_(A.proj_matrices_),
      block_vector_(A.block_vector_, copy_data),
      grid_(A.grid_)
{
    // if(onpe0)cout<<"call ExtendedGridOrbitals(const ExtendedGridOrbitals &A,
    // const bool copy_data)"<<endl;

    assert(A.proj_matrices_ != nullptr);
}

ExtendedGridOrbitals::ExtendedGridOrbitals(const std::string& name,
    const ExtendedGridOrbitals& A, ProjectedMatricesInterface* proj_matrices,
    const bool copy_data)
    : Orbitals(A, copy_data),
      name_(name),
      proj_matrices_(proj_matrices),
      block_vector_(A.block_vector_, copy_data),
      grid_(A.grid_)
{
    assert(proj_matrices != nullptr);

    // setup new projected_matrices object
    proj_matrices_->setup(overlapping_gids_);
}

void ExtendedGridOrbitals::copyDataFrom(const ExtendedGridOrbitals& src)
{
    assert(proj_matrices_ != nullptr);

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
    assert(proj_matrices_ != nullptr);

    if (ct.verbose > 0)
        printWithTimeStamp(
            "ExtendedGridOrbitals::setup()...", (*MPIdata::sout));

    computeGlobalIndexes();

    bool skinny_stencil = !ct.Mehrstellen();

    block_vector_.initialize(overlapping_gids_, skinny_stencil);

    proj_matrices_->setup(overlapping_gids_);

    if (ct.verbose > 0)
        printWithTimeStamp(
            "ExtendedGridOrbitals::setup() done...", (*MPIdata::sout));
}

void ExtendedGridOrbitals::reset(MasksSet* masks, MasksSet* corrmasks,
    std::shared_ptr<LocalizationRegions> lrs)
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

void ExtendedGridOrbitals::initGauss(
    const double rc, const std::shared_ptr<LocalizationRegions> lrs)
{
    assert(numst_ >= 0);
    assert(subdivx_ > 0);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Control& ct     = *(Control::instance());
    if (mmpi.instancePE0() && ct.verbose > 2)
        (*MPIdata::sout) << "Initial orbitals: Gaussians of width " << rc
                         << std::endl;

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
        ORBDTYPE* ipsi               = psi(icolor);
        unsigned int const ipsi_size = numpt_;
        ORBDTYPE* ipsi_host_view     = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(ipsi_size);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            ipsi, ipsi_size, ipsi_host_view);
        MemorySpace::Memory<ORBDTYPE, MemorySpace::Host>::set(
            ipsi_host_view, ipsi_size, 0);

        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const Vector3D& center(lrs->getCenter(icolor));
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
                        const double r = xc.minimage(center, ll, ct.bcWF);
                        if (r < rmax)
                            ipsi_host_view[ix * incx + iy * incy + iz]
                                = static_cast<ORBDTYPE>(exp(-r * r * invrc2));
                        else
                            ipsi_host_view[ix * incx + iy * incy + iz] = 0.;

                        xc[2] += hgrid[2];
                    }
                    xc[1] += hgrid[1];
                }
                xc[0] += hgrid[0];
            }
        }

        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_dev(
            ipsi_host_view, ipsi_size, ipsi);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
            ipsi_host_view);
    }
    resetIterativeIndex();
}

void ExtendedGridOrbitals::initFourier()
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << "Initial orbitals: Fourier " << std::endl;

    const double start0 = grid_.start(0) - grid_.origin(0);
    const double start1 = grid_.start(1) - grid_.origin(1);
    const double start2 = grid_.start(2) - grid_.origin(2);

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
        int index = icolor + 1;
        int kvector[3];
        getkvector(index, cbrtncolors, kvector);

        const double kk[3] = { dk[0] * (double)kvector[0],
            dk[1] * (double)kvector[1], dk[2] * (double)kvector[2] };

        ORBDTYPE* ipsi               = psi(icolor);
        unsigned int const ipsi_size = numpt_;
        ORBDTYPE* ipsi_host_view     = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(ipsi_size);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            ipsi, ipsi_size, ipsi_host_view);
        MemorySpace::Memory<ORBDTYPE, MemorySpace::Host>::set(
            ipsi_host_view, numpt_, 0);

        // TODO this can be done on the GPU with OpenMP
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
                        ipsi_host_view[ix * incx + iy * incy + iz]
                            = 1.
                              - static_cast<ORBDTYPE>(std::cos(kk[0] * x)
                                                      * std::cos(kk[1] * y)
                                                      * std::cos(kk[2] * z));

                        z += hgrid[2];
                    }
                    y += hgrid[1];
                }
                x += hgrid[0];
            }
        }

        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_dev(
            ipsi_host_view, ipsi_size, ipsi);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
            ipsi_host_view);
    }
    resetIterativeIndex();
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

    unsigned int const product_size = numst_ * ldp;
    ORBDTYPE* product_host_view
        = MemorySpace::Memory<ORBDTYPE, memory_space_type>::allocate_host_view(
            product_size);
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
        product, product_size, product_host_view);
    memset(product_host_view, 0, ldp * numst_ * sizeof(ORBDTYPE));

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        unsigned int const phi_size = loc_numpt_ * numst_;
        ORBDTYPE* phi_host_view     = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(phi_size);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            getPsi(0, iloc), phi_size, phi_host_view);

        // TODO this can be done on the GPU
        // Compute product for subdomain iloc
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(loc_numpt_, numst_,
            numst_, 1., phi_host_view, lda_, matrix, numst_, 0.,
            product_host_view + iloc * loc_numpt_, ldp);

        MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
            phi_host_view);
    }
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_dev(
        product_host_view, product_size, product);
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
        product_host_view);

    prod_matrix_tm_.stop();
}

void ExtendedGridOrbitals::multiplyByMatrix(
    const SquareLocalMatrices<MATDTYPE>& matrix, ORBDTYPE* product,
    const int ldp) const
{
    prod_matrix_tm_.start();
    unsigned int const product_size = numst_ * ldp;
    ORBDTYPE* product_host_view
        = MemorySpace::Memory<ORBDTYPE, memory_space_type>::allocate_host_view(
            product_size);
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
        product, product_size, product_host_view);

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        const MATDTYPE* const mat = matrix.getSubMatrix(iloc);

        unsigned int const phi_size = loc_numpt_ * numst_;
        ORBDTYPE* phi_host_view     = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(phi_size);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            getPsi(0, iloc), phi_size, phi_host_view);

        // TODO this can be done on the GPU
        // Compute product for subdomain iloc
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(loc_numpt_, numst_,
            numst_, 1., phi_host_view, lda_, mat, numst_, 0.,
            product_host_view + iloc * loc_numpt_, ldp);
    }

    MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_dev(
        product_host_view, product_size, product);
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
        product_host_view);

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
        const MATDTYPE* const mat   = matrix.getSubMatrix(iloc);
        unsigned int const phi_size = loc_numpt_ * numst_;
        ORBDTYPE* phi_host_view     = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(phi_size);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            getPsi(0, iloc), phi_size, phi_host_view);

        // TODO this can be done on the GPU
        // Compute product for subdomain iloc
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(loc_numpt_, numst_,
            numst_, 1., phi_host_view, lda_, mat, numst_, 0., product,
            loc_numpt_);

        for (int color = 0; color < numst_; color++)
            memcpy(phi_host_view + color * lda_, product + color * loc_numpt_,
                slnumpt);

        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_dev(
            phi_host_view, phi_size, getPsi(0, iloc));
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
            phi_host_view);
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
        unsigned int const phi_size = loc_numpt_ * numst_;
        ORBDTYPE* phi_host_view     = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(phi_size);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            getPsi(0, iloc), phi_size, phi_host_view);

        // TODO this can be done on the GPU
        // Compute loc_numpt_ rows (for subdomain iloc)
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(loc_numpt_, numst_,
            numst_, 1., phi_host_view, lda_, work_matrix, numst_, 0., product,
            loc_numpt_);

        for (int color = 0; color < numst_; color++)
            memcpy(phi_host_view + color * lda_, product + color * loc_numpt_,
                slnumpt);

        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_dev(
            phi_host_view, phi_size, getPsi(0, iloc));
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
            phi_host_view);
    }

    delete[] product;

    prod_matrix_tm_.stop();
}

int ExtendedGridOrbitals::read_hdf5(HDFrestart& h5f_file)
{
    assert(proj_matrices_ != nullptr);

    Control& ct = *(Control::instance());

    hid_t file_id    = h5f_file.file_id();
    std::string name = "Function";
    int ierr         = read_func_hdf5(h5f_file, name);
    if (ierr < 0)
    {
        (*MPIdata::serr)
            << "ExtendedGridOrbitals::read_hdf5(): error in reading " << name
            << ", size=" << name.size() << std::endl;
        return ierr;
    }
    else if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout) << "ExtendedGridOrbitals::read_hdf5(): Read " << ierr
                         << " functions in restart file" << std::endl;
    }

    // Read DM
    if (!ct.fullyOccupied())
    {
        ierr = proj_matrices_->read_dm_hdf5(file_id);
        if (ierr < 0)
        {
            (*MPIdata::serr)
                << "ExtendedGridOrbitals::read_hdf5(): error in reading DM"
                << std::endl;
            return ierr;
        }
    }

    resetIterativeIndex();

    return ierr;
}

int ExtendedGridOrbitals::write_hdf5(
    HDFrestart& h5f_file, const std::string& name)
{
    assert(proj_matrices_ != nullptr);
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

int ExtendedGridOrbitals::write_func_hdf5(
    HDFrestart& h5f_file, const std::string& name)
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
                         << " with precision " << precision << std::endl;
    // loop over global (storage) functions
    for (int color = 0; color < numst_; color++)
    {
        std::string datasetname(getDatasetName(name, color));
        if (onpe0 && ct.verbose > 2)
            (*MPIdata::sout) << "Write " << datasetname << std::endl;

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
                                 << std::endl;
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

        std::vector<int> gids;
        gids.push_back(color);

        if (iwrite)
        {
            writeGids(dset_id, gids);

            // Write the attribute "Lattice parameters" at "Cell origin"
            std::string attname("Lattice parameters");

            // Create the data space for the attribute "Lattice parameters".
            std::vector<double> attr_data(3);
            attr_data[0] = grid_.ll(0);
            attr_data[1] = grid_.ll(1);
            attr_data[2] = grid_.ll(2);

            mgmol_tools::addAttribute2Dataset(
                dset_id, attname.c_str(), attr_data);

            attr_data[0] = grid_.origin(0);
            attr_data[1] = grid_.origin(1);
            attr_data[2] = grid_.origin(2);

            std::string attname2("Cell origin");
            mgmol_tools::addAttribute2Dataset(
                dset_id, attname2.c_str(), attr_data);
        } // iwrite

        unsigned int const psi_size = numpt_;
        ORBDTYPE* psi_host_view     = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(psi_size);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            psi(color), psi_size, psi_host_view);

        int ierr = h5f_file.writeData(
            psi_host_view, filespace, memspace, dset_id, precision);
        if (ierr < 0) return ierr;

        // Close/release resources.
        if (iwrite)
        {
            herr_t status = H5Dclose(dset_id);
            if (status < 0)
            {
                (*MPIdata::serr) << "ExtendedGridOrbitals::write_func_hdf5:"
                                    "H5Dclose failed!!!"
                                 << std::endl;
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
            (*MPIdata::serr) << "H5Sclose filespace failed!!!" << std::endl;
        }
        status = H5Sclose(memspace);
        if (status < 0)
        {
            (*MPIdata::serr) << "H5Sclose memspace failed!!!" << std::endl;
        }
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.barrier();

    return 0;
}

// read all the data sets with names starting with "name"
int ExtendedGridOrbitals::read_func_hdf5(
    HDFrestart& h5f_file, const std::string& name)
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
    hid_t memspace = (h5f_file.active()) ? h5f_file.createMemspace() : 0;

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
                << " PEs" << std::endl;
        }
        else
        {
            (*MPIdata::sout) << "ExtendedGridOrbitals::read_func_hdf5(): Read "
                                "wave functions "
                             << name << " from all tasks..." << std::endl;
        }
    }

    const short precision = ct.restart_info > 3 ? 2 : 1;

    for (int icolor = 0; icolor < numst_; icolor++)
    {
        const std::string datasetname(getDatasetName(name, icolor));

        // check if dataset exists...
        int err_id = h5f_file.dset_exists(datasetname);
        if (h5f_file.gatherDataX()) mmpi.bcast(&err_id, 1);
        if (err_id == 0) break; // dataset does not exists

        if (onpe0 && ct.verbose > 2)
            (*MPIdata::sout) << "Read Dataset " << datasetname
                             << " with precision " << precision << std::endl;

        // Open dataset.
        hid_t dset_id = h5f_file.open_dset(datasetname);
        if (dset_id < 0)
        {
            (*MPIdata::serr)
                << "ExtendedGridOrbitals::read_func_hdf5() --- cannot open "
                << datasetname << std::endl;
            return dset_id;
        }

        herr_t status = h5f_file.readData(buffer, memspace, dset_id, precision);
        if (status < 0)
        {
            (*MPIdata::serr) << "ExtendedGridOrbitals::read_func_hdf5() --- "
                                "H5Dread failed!!!"
                             << std::endl;
            return -1;
        }

        status = h5f_file.close_dset(dset_id);
        if (status < 0)
        {
            return status;
        }

#ifdef HAVE_MAGMA
        ORBDTYPE* buffer_dev
            = MemorySpace::Memory<ORBDTYPE, MemorySpace::Device>::allocate(
                numpt_);
        MemorySpace::copy_to_dev(buffer, numpt_, buffer_dev);
#else
        ORBDTYPE* buffer_dev = buffer;
#endif
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            const int shift = iloc * loc_numpt_;
            block_vector_.assignLocal(icolor, iloc, buffer_dev + shift);
        }
#ifdef HAVE_MAGMA
        MemorySpace::Memory<ORBDTYPE, MemorySpace::Device>::free(buffer_dev);
#endif
    }

    delete[] buffer;
    resetIterativeIndex();

    if (h5f_file.active())
    {
        herr_t status = H5Sclose(memspace);
        if (status < 0)
        {
            (*MPIdata::serr) << "H5Sclose failed!!!" << std::endl;
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

    assert(proj_matrices_ != nullptr);

    matB_tm_.start();
#if DEBUG
    if (onpe0)
        (*MPIdata::sout) << "ExtendedGridOrbitals::computeMatB()" << std::endl;
#endif

    const short bcolor = 32;

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    ORBDTYPE* work = new ORBDTYPE[lda_ * bcolor];
    memset(work, 0, lda_ * bcolor * sizeof(ORBDTYPE));

    ORBDTYPE* const orbitals_psi
        = (numst_ > 0) ? orbitals.block_vector_.vect(0) : nullptr;
    const unsigned int orbitals_psi_size
        = orbitals.block_vector_.get_allocated_size_storage();
    ORBDTYPE* orbitals_psi_host_view
        = MemorySpace::Memory<ORBDTYPE, memory_space_type>::allocate_host_view(
            orbitals_psi_size);
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
        orbitals_psi, orbitals_psi_size, orbitals_psi_host_view);

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
                orbitals_psi_host_view + iloc * loc_numpt_, lda_,
                work + iloc * loc_numpt_, lda_, 0., ssiloc + icolor * numst_,
                numst_);
        }
    }

    MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
        orbitals_psi_host_view);
    delete[] work;

    const double vel = grid_.vel();
    ss.scal(vel);
    proj_matrices_->initializeMatB(ss);

    matB_tm_.stop();
}

// compute <Phi|B|Phi> and its inverse
void ExtendedGridOrbitals::computeBAndInvB(const pb::Lap<ORBDTYPE>& LapOper)
{
    assert(proj_matrices_ != nullptr);

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
        ORBDTYPE* psi        = block_vector_.vect(0);
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            ss.syrk<memory_space_type>(
                iloc, loc_numpt_, psi + iloc * loc_numpt_, lda_);
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
    assert(array != nullptr);
    assert(numst_ != 0);
    assert(grid_.vel() > 0.);
    assert(subdivx_ > 0);

    const ORBDTYPE* const a = transpose ? array : block_vector_.vect(0);
    const ORBDTYPE* const b = transpose ? block_vector_.vect(0) : array;

    const int lda = transpose ? ld : lda_;
    const int ldb = transpose ? lda_ : ld;

    unsigned int const a_size = numpt_ * ss.m();
    ORBDTYPE* a_host_view
        = MemorySpace::Memory<ORBDTYPE, memory_space_type>::allocate_host_view(
            a_size);
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
        const_cast<ORBDTYPE*>(a), a_size, a_host_view);
    unsigned int const b_size = numpt_ * ss.n();
    ORBDTYPE* b_host_view
        = MemorySpace::Memory<ORBDTYPE, memory_space_type>::allocate_host_view(
            b_size);
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
        const_cast<ORBDTYPE*>(b), b_size, b_host_view);

#ifdef USE_MP
    // use temporary float data for matrix ss
    LocalMatrices<ORBDTYPE> ssf(ss.subdiv(), ss.m(), ss.n());
#else
    LocalMatrices<ORBDTYPE>& ssf(ss);
#endif
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        ssf.gemm(iloc, loc_numpt_, a_host_view + iloc * loc_numpt_, lda,
            b_host_view + iloc * loc_numpt_, ldb);
    }
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
        a_host_view);
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
        b_host_view);
#ifdef USE_MP
    ss.copy(ssf);
#endif

    ss.scal(grid_.vel());
}

void ExtendedGridOrbitals::computeDiagonalElementsDotProduct(
    const ExtendedGridOrbitals& orbitals, std::vector<DISTMATDTYPE>& ss) const
{
    assert(numst_ > 0);
    assert(grid_.vel() > 0.);

    for (int icolor = 0; icolor < numst_; icolor++)
    {
        ss[icolor] = 0.;
        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            double alpha
                = LinearAlgebraUtils<memory_space_type>::MPdot(loc_numpt_,
                    orbitals.getPsi(icolor, iloc), getPsi(icolor, iloc));

            ss[icolor] += (DISTMATDTYPE)(alpha * grid_.vel());
        }
    }
    std::vector<DISTMATDTYPE> tmp(ss);
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&tmp[0], &ss[0], numst_, MPI_SUM);
}

void ExtendedGridOrbitals::computeGram(
    dist_matrix::DistMatrix<DISTMATDTYPE>& gram_mat)
{
    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    getLocalOverlap(ss);

    LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();

    gram_mat.clear();

    sl2dm->accumulate(ss, gram_mat);
}

void ExtendedGridOrbitals::computeGram(const ExtendedGridOrbitals& orbitals,
    dist_matrix::DistMatrix<DISTMATDTYPE>& gram_mat)
{
    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    getLocalOverlap(orbitals, ss);

    // make a DistMatrix out of ss
    LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();

    gram_mat.clear();

    sl2dm->accumulate(ss, gram_mat);
}

// compute the lower-triangular part of the overlap matrix
void ExtendedGridOrbitals::computeGram(const int verbosity)
{
    assert(proj_matrices_ != nullptr);

    overlap_tm_.start();

#ifdef PRINT_OPERATIONS
    if (onpe0)
        (*MPIdata::sout) << "ExtendedGridOrbitals::computeGram()" << std::endl;
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
    assert(proj_matrices_ != nullptr);

    computeGram(verbosity);

    /* Compute inverse of Gram matrix */
    proj_matrices_->computeInvS();
}

void ExtendedGridOrbitals::checkCond(const double tol, const bool flag_stop)
{
    assert(proj_matrices_ != nullptr);

    proj_matrices_->checkCond(tol, flag_stop);
}

double ExtendedGridOrbitals::dotProductWithDM(
    const ExtendedGridOrbitals& orbitals)
{
    assert(proj_matrices_ != nullptr);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    computeLocalProduct(orbitals, ss);

    return proj_matrices_->dotProductWithDM(ss);
}

double ExtendedGridOrbitals::dotProductWithInvS(
    const ExtendedGridOrbitals& orbitals)
{
    assert(proj_matrices_ != nullptr);

    SquareLocalMatrices<MATDTYPE> ss(subdivx_, numst_);

    computeLocalProduct(orbitals, ss);

    return proj_matrices_->dotProductWithInvS(ss);
}

double ExtendedGridOrbitals::dotProductDiagonal(
    const ExtendedGridOrbitals& orbitals)
{
    assert(proj_matrices_ != nullptr);

    std::vector<DISTMATDTYPE> ss(numst_);
    computeDiagonalElementsDotProduct(orbitals, ss);
    return proj_matrices_->getTraceDiagProductWithInvS(ss);
}

double ExtendedGridOrbitals::dotProductSimple(
    const ExtendedGridOrbitals& orbitals)
{
    assert(proj_matrices_ != nullptr);

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
                         << std::endl;
        Control& ct = *(Control::instance());
        ct.global_exit(2);
    }

    dot_product_tm_.stop();

    return dot;
}

dist_matrix::DistMatrix<DISTMATDTYPE> ExtendedGridOrbitals::product(
    const ExtendedGridOrbitals& orbitals, const bool transpose)
{
    assert(numst_ > 0);
    assert(subdivx_ > 0);
    assert(subdivx_ < 1000);

    return product(orbitals.psi(0), numst_, orbitals.lda_, transpose);
}

dist_matrix::DistMatrix<DISTMATDTYPE> ExtendedGridOrbitals::product(
    const ORBDTYPE* const array, const int ncol, const int lda,
    const bool transpose)
{
    assert(lda > 1);

    dot_product_tm_.start();

    LocalMatrices<MATDTYPE> ss(subdivx_, numst_, ncol);

    computeLocalProduct(array, lda, ss, transpose);

    LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();

    dist_matrix::DistMatrix<DISTMATDTYPE> tmp("tmp", numst_, numst_);
    sl2dm->accumulate(ss, tmp);

    dot_product_tm_.stop();

    return tmp;
}

void ExtendedGridOrbitals::orthonormalizeLoewdin(const bool overlap_uptodate,
    SquareLocalMatrices<MATDTYPE>* matrixTransform, const bool update_matrices)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "ExtendedGridOrbitals::orthonormalizeLoewdin()"
                         << std::endl;

    if (!overlap_uptodate) computeGram(0);

    SquareLocalMatrices<MATDTYPE>* localP = matrixTransform;
    if (matrixTransform == nullptr)
        localP = new SquareLocalMatrices<MATDTYPE>(subdivx_, numst_);

    incrementIterativeIndex();

    bool multbymat = false;
#ifdef HAVE_MAGMA
    // try with ReplicatedMatrix first
    {
        ProjectedMatrices<ReplicatedMatrix>* projmatrices
            = dynamic_cast<ProjectedMatrices<ReplicatedMatrix>*>(
                proj_matrices_);
        if (projmatrices)
        {
            projmatrices->computeLoewdinTransform(
                *localP, getIterativeIndex(), update_matrices);
            multiplyByMatrix(*localP);

            projmatrices->setGram2Id(getIterativeIndex());

            multbymat = true;
        }
    }
#endif
    if (!multbymat)
    {
        ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>* projmatrices
            = dynamic_cast<
                ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>*>(
                proj_matrices_);
        if (projmatrices)
        {
            projmatrices->computeLoewdinTransform(
                *localP, getIterativeIndex(), update_matrices);
            multiplyByMatrix(*localP);

            projmatrices->setGram2Id(getIterativeIndex());
        }
    }

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
                         << st1 << " and " << st2 << std::endl;
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
                         << "," << overlap[2] << std::endl;
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
    std::vector<std::vector<double>>& inv_norms2) const
{
    std::vector<double> diagS(numst_);

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
    std::vector<double> diagS(numst_);

    computeDiagonalElementsDotProduct(*this, diagS);

    for (int color = 0; color < numst_; color++)
    {
#ifdef DEBUG
        if (onpe0 && ct.verbose > 2)
            for (int i = 0; i < numst_; i++)
                (*MPIdata::sout)
                    << "i=" << i << ", diagS[i]=" << diagS[i] << std::endl;
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
    assert(lda_ >= loc_numpt_);

    SquareLocalMatrices<MATDTYPE> lmatrix(subdivx_, numst_);

    if (numst_ != 0) computeLocalProduct(array, lda, lmatrix, false);

#ifdef DEBUG
    (*MPIdata::sout) << "ExtendedGridOrbitals::projectOut()" << std::endl;
    (*MPIdata::sout) << "Product before projection" << std::endl;
    pmatrix.print((*MPIdata::sout));
#endif
    proj_matrices_->applyInvS(lmatrix);

    ORBDTYPE* tproduct = new ORBDTYPE[loc_numpt_ * numst_];
    memset(tproduct, 0, loc_numpt_ * numst_ * sizeof(ORBDTYPE));

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        unsigned int const phi_size = loc_numpt_ * numst_;
        ORBDTYPE* phi_host_view     = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(phi_size);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            getPsi(0, iloc), phi_size, phi_host_view);

        MATDTYPE* localMat_iloc = lmatrix.getSubMatrix(iloc);

        // TODO this can be done on the GPU
        // Compute loc_numpt_ rows (for subdomain iloc)
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(loc_numpt_, numst_,
            numst_, 1., phi_host_view, lda_, localMat_iloc, numst_, 0.,
            tproduct, loc_numpt_);

        MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
            phi_host_view);

        ORBDTYPE* parray               = array + iloc * loc_numpt_;
        unsigned int const parray_size = numst_ * lda;
        ORBDTYPE* parray_host_view     = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(parray_size);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            parray, parray_size, parray_host_view);

        double minus = -1. * scale;
        for (int j = 0; j < numst_; j++)
            LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(loc_numpt_, minus,
                tproduct + j * loc_numpt_, parray_host_view + j * lda);

        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_dev(
            parray_host_view, parray_size, parray);

        MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
            parray_host_view);
    }

    delete[] tproduct;
}

void ExtendedGridOrbitals::initRand()
{
    Control& ct = *(Control::instance());

    const unsigned dim[3] = { grid_.dim(0), grid_.dim(1), grid_.dim(2) };

    std::vector<double> xrand(grid_.gdim(0));
    std::vector<double> yrand(grid_.gdim(1));
    std::vector<double> zrand(grid_.gdim(2));

    const int loc_length = dim[0] / subdivx_;
    assert(loc_length > 0);
    assert(static_cast<unsigned int>(loc_length) <= dim[0]);

    const int xoff = grid_.istart(0);
    const int yoff = grid_.istart(1);
    const int zoff = grid_.istart(2);

    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << " Initialize " << numst_
                         << " random global functions" << std::endl;

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

        unsigned int const size  = loc_numpt_;
        ORBDTYPE* psi_state_view = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(size);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            psi(istate), size, psi_state_view);

        for (short iloc = 0; iloc < subdivx_; iloc++)
        {
            for (int ix = loc_length * iloc; ix < loc_length * (iloc + 1); ix++)
                for (unsigned int iy = 0; iy < dim[1]; iy++)
                    for (unsigned int iz = 0; iz < dim[2]; iz++)
                    {
                        const double alpha = xrand[xoff + ix] * yrand[yoff + iy]
                                             * zrand[zoff + iz];

                        psi_state_view[ix * incx + iy * incy + iz]
                            = alpha * alpha;

                        assert((ix * incx + iy * incy + iz)
                               < static_cast<unsigned int>(lda_));
                    }
        }
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_dev(
            psi_state_view, size, psi(istate));
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
            psi_state_view);
    }

    resetIterativeIndex();
}

template<>
void ExtendedGridOrbitals::addDotWithNcol2Matrix(
    ExtendedGridOrbitals& Apsi, dist_matrix::DistMatrix<double>& matrix) const
{
    addDot_tm_.start();

    assert(numst_ > 0);

    const double vel = grid_.vel();

    // replicated matrix
    const int size_work = numst_ * numst_;
    std::vector<double> work(size_work);
    memset(work.data(), 0, size_work * sizeof(double));

    unsigned int const block_vector_size = numpt_ * numst_;
    ORBDTYPE* block_vector_host_view
        = MemorySpace::Memory<ORBDTYPE, memory_space_type>::allocate_host_view(
            block_vector_size);
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
        block_vector_.vect(0), block_vector_size, block_vector_host_view);

    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        unsigned int const phi_size = loc_numpt_ * numst_;
        ORBDTYPE* phi_host_view     = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(phi_size);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            Apsi.getPsi(0, iloc), phi_size, phi_host_view);

        // TODO this can be done on the GPU
        MPgemmTN(numst_, numst_, loc_numpt_, vel,
            block_vector_host_view + iloc * loc_numpt_, lda_, phi_host_view,
            lda_, 1., work.data(), numst_);

        MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
            phi_host_view);
    }
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
        block_vector_host_view);

    std::vector<double> work2(size_work);
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(work.data(), work2.data(), numst_ * numst_, MPI_SUM);

    // replicated -> DistMatrix
    matrix.add(work2.data(), numst_);

    addDot_tm_.stop();
}

#ifdef HAVE_MAGMA
template<>
void ExtendedGridOrbitals::addDotWithNcol2Matrix(
    ExtendedGridOrbitals& Apsi, ReplicatedMatrix& matrix) const
{
    addDot_tm_.start();

    magma_trans_t magma_transa = magma_trans_const('t');
    magma_trans_t magma_transb = magma_trans_const('n');

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magmablas_dgemm(magma_transa, magma_transb, numst_, numst_, numpt_, 1.,
        block_vector_.vect(0), lda_, Apsi.getPsi(0), lda_, 0.,
        matrix.data(), matrix.ld(), magma_singleton.queue_);

    addDot_tm_.stop();
}
#endif

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

void ExtendedGridOrbitals::printTimers(std::ostream& os)
{
    matB_tm_.print(os);
    invBmat_tm_.print(os);
    overlap_tm_.print(os);
    dot_product_tm_.print(os);
    addDot_tm_.print(os);
    prod_matrix_tm_.print(os);
    assign_tm_.print(os);
    normalize_tm_.print(os);
    axpy_tm_.print(os);
}

void ExtendedGridOrbitals::initWF(
    const std::shared_ptr<LocalizationRegions> lrs)
{
    Control& ct = *(Control::instance());

    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << " Initialize wave functions ..." << std::endl;
    }
    switch (ct.init_type)
    {
        case 1:
            if (onpe0 && ct.verbose > 1)
            {
                (*MPIdata::sout) << " with Gaussian functions..." << std::endl;
            }
            initGauss(ct.init_rc, lrs);
            break;
        case 2:
            if (onpe0 && ct.verbose > 1)
            {
                (*MPIdata::sout) << " with Fourier basis ..." << std::endl;
            }
            initFourier();
            break;
        default:
            if (onpe0 && ct.verbose > 2)
            {
                (*MPIdata::sout) << " with random values ..." << std::endl;
            }
            initRand();

            if (ct.globalColoring())
            {
                // smooth out random functions
                pb::Laph4M<ORBDTYPE> myoper(grid_);
                pb::GridFunc<ORBDTYPE> gf_work(
                    grid_, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
                pb::GridFunc<ORBDTYPE> gf_psi(
                    grid_, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);

                if (onpe0 && ct.verbose > 2)
                    (*MPIdata::sout)
                        << " Apply B to initial wave functions" << std::endl;
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
            << " Normalize or Orthonormalize initial wave functions"
            << std::endl;
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
        (*MPIdata::sout) << "ExtendedGridOrbitals::init_wf() done" << std::endl;
#endif
}

template void ExtendedGridOrbitals::setDataWithGhosts(
    pb::GridFuncVector<float, memory_space_type>* data_wghosts);
template void ExtendedGridOrbitals::setDataWithGhosts(
    pb::GridFuncVector<double, memory_space_type>* data_wghosts);

template void ExtendedGridOrbitals::setPsi(
    const pb::GridFunc<float>& gf_work, const int ist);
template void ExtendedGridOrbitals::setPsi(
    const pb::GridFunc<double>& gf_work, const int ist);

template void ExtendedGridOrbitals::setPsi(
    const pb::GridFuncVector<float, memory_space_type>& gf_work);
template void ExtendedGridOrbitals::setPsi(
    const pb::GridFuncVector<double, memory_space_type>& gf_work);
