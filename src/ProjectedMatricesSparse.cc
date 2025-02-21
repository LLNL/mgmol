// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ProjectedMatricesSparse.h"
#include <cassert>
#include <fstream>
#include <iomanip>

#include "Control.h"
#include "HDFrestart.h"
#include "MGmol_MPI.h"
#include "VariableSizeMatrix.h"
#include "random.h"
#include "tools.h"

#define ORBITAL_OCCUPATION 2.0
#define RY2EV 13.605804

Timer ProjectedMatricesSparse::compute_inverse_tm_(
    "ProjectedMatricesSparse::computeInverse");
Timer ProjectedMatricesSparse::compute_invB_tm_(
    "ProjectedMatricesSparse::computeInvB");
Timer ProjectedMatricesSparse::init_gram_matrix_tm_(
    "ProjectedMatricesSparse::initialize_Gram_Matrix");
Timer ProjectedMatricesSparse::update_theta_tm_(
    "ProjectedMatricesSparse::updateTheta");
Timer ProjectedMatricesSparse::update_submatT_tm_(
    "ProjectedMatricesSparse::updateSubmatT");
Timer ProjectedMatricesSparse::update_submatX_tm_(
    "ProjectedMatricesSparse::updateSubmatX");
Timer ProjectedMatricesSparse::eigsum_tm_("ProjectedMatricesSparse::eigsum");
Timer ProjectedMatricesSparse::consolidate_H_tm_(
    "ProjectedMatricesSparse::consolidate_sH");
Timer ProjectedMatricesSparse::eig_interval_tm_(
    "ProjectedMatrices::computeEigenInterval");

ProjectedMatricesSparse::ProjectedMatricesSparse(const int ndim,
    const double width, std::shared_ptr<LocalizationRegions> lrs,
    ClusterOrbitals* local_cluster)
    : ProjectedMatricesInterface(false, width)
{
    assert(lrs);

    dim_     = ndim;
    min_val_ = 0.25;

    lrs_           = lrs;
    local_cluster_ = local_cluster;

    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    distributor_invS_
        = new DataDistribution("invS", lrs_->max_radii(), myPEenv, domain);

    invS_ = nullptr;

    dm_ = nullptr;

    sH_      = nullptr;
    matHB_   = nullptr;
    submatT_ = nullptr;
    localX_  = nullptr;
    localT_  = nullptr;

    isDataSetup_ = false;
    /* Initialize eigenvalue0_ = 0. for use in OrbitalsPreconditioning.
     * This is not set in the sparse case and perhaps
     * could be approximated with the diagonal of the
     * Gram matrix. It could also be removed if not
     * needed for computing gamma in OrbitalsPreconditioning.
     */
    eigenvalue0_ = 0.;
}

void ProjectedMatricesSparse::clearData()
{
    //    if(isDataSetup_)
    //    {
    //       delete invS_; invS_=0;
    //       delete dm_; dm_ = 0;
    delete localX_;
    localX_ = nullptr;
    delete localT_;
    localT_ = nullptr;
    delete sH_;
    sH_ = nullptr;
    delete matHB_;
    matHB_ = nullptr;
    delete submatT_;
    submatT_ = nullptr;

    isDataSetup_ = false;
    //    }
}

ProjectedMatricesSparse::~ProjectedMatricesSparse()
{
    assert(invS_ != nullptr);
    assert(dm_ != nullptr);
    assert(localX_ != nullptr);
    assert(localT_ != nullptr);
    assert(submatT_ != nullptr);
    assert(sH_ != nullptr);
    assert(matHB_ != nullptr);
    assert(distributor_invS_ != nullptr);

    delete distributor_invS_;
    distributor_invS_ = nullptr;
    delete invS_;
    invS_ = nullptr;
    delete dm_;
    dm_ = nullptr;
    clearData();
}

void ProjectedMatricesSparse::setup(
    const std::vector<std::vector<int>>& global_indexes)
{
    // assert( (short)global_indexes.size()>0 );
    // assert( (short)global_indexes[0].size()>0 );
    Control& ct = *(Control::instance());

    setupBase(global_indexes.size(), global_indexes[0].size());

    global_indexes_ = global_indexes;

    locvars_.clear();

    // clear old data
    clearData();
    // clear invS and DM data if MD
    // Ideally MD would be calling a different setup routine
    // from LocGridOrbitals, that resets all data
    if (ct.AtomsDynamic() == AtomsDynamicType::MD)
    {
        delete invS_;
        invS_ = nullptr;
        delete dm_;
        dm_ = nullptr;
    }

    // printf("setup is called ...\n");
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();
    int myid;
    MPI_Comm_rank(comm, &myid);

    std::vector<int> locfcns;
    (*lrs_).getLocalSubdomainIndices(locfcns);

    if (dim_ > 0)
    {
        /* setup matrix data */
        locvars_ = (*lrs_).getOverlapGids();

        lsize_ = locvars_.size();

        // reset invS and DM data
        if (invS_ == nullptr || dm_ == nullptr
            || ct.AtomsDynamic() == AtomsDynamicType::MD)
        {
            invS_ = new ShortSightedInverse(lrs_, locvars_, local_cluster_);
            dm_   = new DensityMatrixSparse(lrs_, dim_, locvars_);
        }

        matHB_ = new VariableSizeMatrix<sparserow>("HB", lsize_);
        // estimate size of table needed for efficient access to elements of sH
        sH_ = new VariableSizeMatrix<sparserow>("sH", 4096);
        /* initialize Sparse H matrix -- this is necessary for efficient data
         * distribution */
        (*sH_).setupSparseRows(locvars_);

        submatT_ = new VariableSizeMatrix<sparserow>("Theta", lsize_);

        localX_ = new SquareLocalMatrices<MATDTYPE, memory_space_type>(
            subdiv_, chromatic_number_);
        localT_ = new SquareLocalMatrices<MATDTYPE, MemorySpace::Host>(
            subdiv_, chromatic_number_);
    }
    // get min and max of functions centered on local processors
    ///*
    int locmin, locmax;
    locmin = locfcns.size();
    locmax = locmin;
    mmpi.allreduce(&locmin, 1, MPI_MIN);
    mmpi.allreduce(&locmax, 1, MPI_MAX);

    if (mmpi.instancePE0() && ct.verbose > 0)
    {
        std::cout << "Max. number of locally centered functions: " << locmax
                  << std::endl;
        std::cout << "Min. number of locally centered functions: " << locmin
                  << std::endl;
    }
    //*/

    isDataSetup_ = true;
}

void ProjectedMatricesSparse::updateSubMatT()
{
    update_submatT_tm_.start();

    updateLocalMat(*submatT_, localT_);

    update_submatT_tm_.stop();
}

/* compute theta = invB * Hij */
void ProjectedMatricesSparse::updateTheta()
{
    assert(invS_ != nullptr);

    update_theta_tm_.start();

    (*submatT_).setupSparseRows(locvars_);
    /* get local contribution to theta */
    (*invS_).invSmultB(sH_, (*submatT_));

    /* gather over neighboring processes */
    double spread_radius     = (*lrs_).max_radii();
    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    DataDistribution distributor("Theta", spread_radius, myPEenv, domain);

    distributor.updateLocalRows((*submatT_));

    update_theta_tm_.stop();
}

double ProjectedMatricesSparse::getExpectationH()
{
    assert(invS_ != nullptr);
    assert(matHB_ != nullptr);
    //   assert(matHB_->n() != 0);

    return dm_->getTraceDotProductWithMat(matHB_);
}

double ProjectedMatricesSparse::getExpectationMat(
    VariableSizeMatrix<sparserow>* mat)
{
    assert(invS_ != nullptr);
    assert(mat != nullptr);
    assert(mat->n() != 0);

    return dm_->getTraceDotProductWithMat(mat);
}

void ProjectedMatricesSparse::consolidateOrbitalsOverlapMat(
    VariableSizeMatrix<sparserow>& mat)
{
    consolidate_H_tm_.start();

    std::vector<int> locfcns;
    lrs_->getLocalSubdomainIndices(locfcns);

    // update locally centered row data
    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    // consolidate locally centered data
    DataDistribution distributor1(
        "overlap", lrs_->max_radii(), myPEenv, domain);
    distributor1.augmentLocalData(mat, false);

    // now gather updated data from neighbors.
    // gather from neighbors that contain functions that overlap with
    // locally centered functions
    Control& ct(*Control::instance());
    DataDistribution distributor("overlap2", ct.spread_radius, myPEenv, domain);

    //    mat.consolidate(locfcns, distributor);
    distributor.consolidateMatrix(locfcns, mat);

    consolidate_H_tm_.stop();
}

// Gather data to initialize matH and matHB
void ProjectedMatricesSparse::consolidateH()
{
    // matrix size is at least as large as the number of overlapping orbitals
    // can be larger with nonlocal projectors that "couple"
    // non-overlapping orbitals
    assert(sH_->n() >= (int)locvars_.size());

    consolidate_H_tm_.start();

    std::vector<int> locfcns;
    lrs_->getLocalSubdomainIndices(locfcns);

    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };

    /* gather data for matH amd matHB */
    DataDistribution distributor_sH(
        "prodH", 3. * lrs_->max_radii(), myPEenv, domain);
    distributor_sH.augmentLocalData((*sH_), false);

    // now gather updated data from neighbors.
    // gather from neighbors that contain functions that overlap with locally
    // centered functions
    Control& ct(*(Control::instance()));
    DataDistribution distributor("H", ct.spread_radius, myPEenv, domain);

    distributor.consolidateMatrix(locfcns, (*sH_));

    matHB_->copyData((*sH_), lsize_);

    consolidate_H_tm_.stop();
}

// computes the trace of the dotproduct with invS
// assumes matrix ss contains partial contributions
double ProjectedMatricesSparse::dotProductWithInvS(
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& local_product)
{
    return computeTraceInvSmultMat(local_product, true);
}

// computes the trace of the dotproduct with DM
// assumes matrix ss contains partial contributions
double ProjectedMatricesSparse::dotProductWithDM(
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& local_product)
{
    double trace = computeTraceDMmultMat(local_product, true);
    // scale by occupation number
    // NOTE: Ideally this scaling should be done outside this function
    MGmol_MPI& mmpi           = *(MGmol_MPI::instance());
    double orbital_occupation = mmpi.nspin() > 1 ? 1. : 2.;
    trace /= orbital_occupation;

    return trace;
}

double ProjectedMatricesSparse::dotProductSimple(
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& /*local_product*/)
{
    std::cerr << "ERROR: ProjectedMatricesSparse::dotProductSimple() \
               not implemented!!!"
              << std::endl;

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Abort(mmpi.commSameSpin(), EXIT_FAILURE);

    return -1.;
}

void ProjectedMatricesSparse::printMatrices(std::ostream& os) const
{
    printS(os);
    printH(os);
    printTheta(os);
    printHB(os);
}

double ProjectedMatricesSparse::getNel() const { return 2. * dim_; }

double ProjectedMatricesSparse::getEigSum()
{
    if (dim_ == 0) return 0.;
    assert(invS_ != nullptr);
    eigsum_tm_.start();
    /* local functions for computing trace */
    std::vector<int> locfcns;
    (*lrs_).getLocalSubdomainIndices(locfcns);

    /* compute trace */
    double trace = 0.0;
    for (std::vector<int>::iterator itr = locfcns.begin(); itr != locfcns.end();
         ++itr)
    {
        trace += (*invS_).invSmultB_ij(matHB_, *itr, *itr);
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();

    double val = 0.0;
    MPI_Allreduce(&trace, &val, 1, MPI_DOUBLE, MPI_SUM, comm);
    eigsum_tm_.stop();
    return val;
}

double ProjectedMatricesSparse::getLinDependent2states(
    int& st1, int& st2, const bool print_flag) const
{
    assert((*invS_).gramMat() != nullptr);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();

    /* get pointer to new Gram matrix */
    VariableSizeMatrix<sparserowtab>* mat = (*invS_).gramMat();

    /* get size of submatrix and set its extents */
    const int bsize  = locvars_.size();
    const int bnzmax = (*mat).getNzmaxSubmat(0, bsize - 1);

    /* construct Linear Solver matrix */
    LinearSolverMatrix<lsdatatype> LSMat(bsize, bnzmax);
    LSMat.initSquareMat((*mat), false);

    /* get smallest eigenvalue of LSMat, return corresponding eigenvector */
    LinearSolver solver;
    std::vector<double> evec; //(bsize, 0.);
    double eigmin = solver.computeEigMin(LSMat, evec);
    /* get linearly dependent states */
    std::vector<double>::iterator idx1
        = std::max_element(evec.begin(), evec.end());
    std::vector<double>::iterator idx2
        = std::min_element(evec.begin(), evec.end());
    std::vector<double>::iterator idx
        = (fabs(*idx1) > fabs(*idx2)) ? idx1 : idx2;

    assert(idx != evec.end());
    const int locst1 = distance(evec.begin(), idx);
    assert(locst1 >= 0);

#if 0
   vector<double> evec_save(evec);
   //set largest entry to zero and search for second largest entry
   *idx = 0.;
   idx1 = std::max_element(evec.begin(), evec.end());
   idx2 = std::min_element(evec.begin(), evec.end());
   idx = (fabs(*idx1) > fabs(*idx2)) ? idx1 : idx2;
   assert( idx!=evec.end() );
   const int locst2 = distance( evec.begin(), idx );
   assert( locst2>=0 );
   const double vmax2=*idx;
   int states[2] = {locvars_[locst1], locvars_[locst2]};
   evec=evec_save;
#else
    // get second index by looking for largest off-diagonal element in row
    // locvars_[locst1]
    double vmax2  = 0.;
    int states[2] = { locvars_[locst1], -1 };
    // const int locst2 = mat->getMaxAbsOffDiagonalRowEntry(locst1,vmax2);
    states[1] = mat->getMaxAbsOffDiagonalRowEntry(locvars_[locst1], vmax2);
#endif

    // get global minimum
    int myid;
    MPI_Comm_rank(comm, &myid);
    typeDoubleInt in;
    typeDoubleInt out;
    in.val  = eigmin;
    in.rank = myid;
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);

    // get global states and broadcast st1 and st2 from PE with min eigmin
    MPI_Bcast(states, 2, MPI_INT, out.rank, comm);

    st1 = states[0];
    st2 = states[1];

    if (myid == out.rank && print_flag)
    {
        std::cout << "" << std::endl;
        std::cout
            << "ProjectedMatricesSparse::getLinDependent2states(): MPI task "
            << myid << std::endl;
        std::cout << "eigenvector size: " << evec.size() << std::endl;
        // std::cout<<"locst1="<<locst1<<", locst2="<<locst2<<endl;
        std::cout << "off-diagonal element: " << vmax2 << std::endl;
        std::cout << "2x2 submatrix associated with small eigenvalue:"
                  << std::endl;
        mat->printMatBlock2(st1, st2, std::cout);
        std::vector<double> tmp(bsize);
        LSMat.matvec(&evec[0], &tmp[0]);
        int one      = 1;
        double norme = DNRM2(&bsize, &evec[0], &one);
        double norma = DNRM2(&bsize, &tmp[0], &one);
        std::cout << "eigenvalue: " << out.val << std::endl;
        std::cout << "Norm S*v: " << norma / norme << std::endl;
    }
    if (print_flag) mmpi.barrier();

    assert(st1 >= 0);
    assert(st2 >= 0);

    return (out.val);
}

void ProjectedMatricesSparse::printGramMatrix2states(
    const int st1, const int st2, std::ostream& os) const
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();

    int myid;
    MPI_Comm_rank(comm, &myid);
    typeDoubleInt in;
    typeDoubleInt out;
    in.val = 0.;
    if (find(locvars_.begin(), locvars_.end(), st1) != locvars_.end())
        in.val += 1.;
    if (find(locvars_.begin(), locvars_.end(), st2) != locvars_.end())
        in.val += 1.;
    in.rank = myid;
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);

    if (myid == out.rank)
    {
        VariableSizeMatrix<sparserowtab>* mat = invS_->gramMat();
        mat->printMatBlock2(st1, st2, os);
    }
}

void ProjectedMatricesSparse::printTimers(std::ostream& os)
{
    compute_inverse_tm_.print(os);
    compute_invB_tm_.print(os);
    update_theta_tm_.print(os);
    update_submatX_tm_.print(os);
    update_submatT_tm_.print(os);
    init_gram_matrix_tm_.print(os);
    eigsum_tm_.print(os);
    consolidate_H_tm_.print(os);
}

void ProjectedMatricesSparse::updateLocalMat(
    const VariableSizeMatrix<sparserow>& submatM,
    SquareLocalMatrices<MATDTYPE, MemorySpace::Host>* localM)
{
    // return if matrix size is zero
    if (localM->n() == 0) return;

    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        MATDTYPE* localM_iloc = localM->getRawPtr(iloc);
        for (int icolor = 0; icolor < chromatic_number_; icolor++)
        {
            const int st1 = global_indexes_[iloc][icolor];
            if (st1 == -1) continue;
            for (int jcolor = 0; jcolor < chromatic_number_; jcolor++)
            {
                const int st2 = global_indexes_[iloc][jcolor];
                if (st2 == -1) continue;
                localM_iloc[icolor + chromatic_number_ * jcolor]
                    = submatM.get_value(st1, st2);
            }
        }
    }
}

double ProjectedMatricesSparse::computeTraceInvSmultMat(
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& mat,
    const bool consolidate)
{
    Control& ct = *(Control::instance());
    VariableSizeMatrix<sparserow> vsmat("VS", lsize_);
    vsmat.setupSparseRows(locvars_);
    vsmat.insertMatrixElements(mat, global_indexes_, ct.numst);

    if (consolidate)
    {
        // consolidate only locally centered data since trace is needed
        Mesh* mymesh             = Mesh::instance();
        const pb::Grid& mygrid   = mymesh->grid();
        const pb::PEenv& myPEenv = mymesh->peenv();
        double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
        DataDistribution distributor(
            "traceSqLo", (*lrs_).max_radii(), myPEenv, domain);
        distributor.augmentLocalData(vsmat, false);
    }

    return invS_->getTraceDotProductWithInvS(&vsmat);
}

double ProjectedMatricesSparse::computeTraceDMmultMat(
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& mat,
    const bool consolidate)
{
    Control& ct = *(Control::instance());
    VariableSizeMatrix<sparserow> vsmat("VS", lsize_);
    vsmat.setupSparseRows(locvars_);
    vsmat.insertMatrixElements(mat, global_indexes_, ct.numst);

    if (consolidate)
    {
        // consolidate only locally centered data since trace is needed
        Mesh* mymesh             = Mesh::instance();
        const pb::Grid& mygrid   = mymesh->grid();
        const pb::PEenv& myPEenv = mymesh->peenv();
        double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
        DataDistribution distributor(
            "traceSqLo", (*lrs_).max_radii(), myPEenv, domain);
        distributor.augmentLocalData(vsmat, false);
    }

    return dm_->getTraceDotProductWithMat(&vsmat);
}

double ProjectedMatricesSparse::computeTraceInvSmultMat(
    VariableSizeMatrix<sparserow>& vsmat, const bool consolidate)
{
    if (consolidate)
    {
        // consolidate only locally centered data since trace is needed
        Mesh* mymesh             = Mesh::instance();
        const pb::Grid& mygrid   = mymesh->grid();
        const pb::PEenv& myPEenv = mymesh->peenv();
        double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
        DataDistribution distributor(
            "trace", (*lrs_).max_radii(), myPEenv, domain);
        distributor.augmentLocalData(vsmat, false);
    }

    return invS_->getTraceDotProductWithInvS(&vsmat);
}

void ProjectedMatricesSparse::applyInvS(
    SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& mat)
{
    Control& ct = *(Control::instance());
    mat.transpose();

    // convert mat into variablesizematrix object
    VariableSizeMatrix<sparserow> vsmat("vs", lsize_);
    vsmat.setupSparseRows(locvars_);
    vsmat.insertMatrixElements(mat, global_indexes_, ct.numst);

    // consolidate mat to gather data for matmult
    consolidateOrbitalsOverlapMat(vsmat);

    // compute matmult with invS to get locally centered rows of product
    VariableSizeMatrix<sparserow> pmat("p", lsize_);
    pmat.setupSparseRows(locvars_);
    invS_->invSmultB(&vsmat, pmat, false);

    /* gather over neighboring processes */
    double spread_radius     = lrs_->max_radii();
    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    DataDistribution distributor("applyInvS", spread_radius, myPEenv, domain);

    distributor.updateLocalRows(pmat);

    // gather result into square local matrix
    mat.reset();
    updateLocalMat(pmat, &mat);
}

//-----------------------------------------------
#include "linear_algebra/lapack_c.h"

/* Use the power method to compute the extents of the spectrum of the
 * generalized eigenproblem. In order to use a residual-based convergence
 * criterion in an efficient way, we delay normalization of the vectors to avoid
 * multiple matvecs. NOTE: We are only interested in the eigenvalues, so the
 * final eigenvector may not be normalized.
 */

void ProjectedMatricesSparse::computeGenEigenInterval(
    std::vector<double>& interval, const int maxits, const double pad)
{
    srand(13579);
    assert(invS_ != nullptr);
    assert(sH_->n() == invS_->gramMat()->n());

    eig_interval_tm_.start();

    interval.clear();

    VariableSizeMatrix<sparserow> mat(*sH_);

    // use the power method to get the eigenvalue interval
    const int m      = mat.n(); // number of (local) rows
    const double one = 1., zero = 0.;

    // define static variables and shift for matrices to speed up convergence
    static double shft = 0.; // shft is initially zero
    // shift
    mat.axpy(shft, *(invS_->gramMat()));

    // initialize random vectors for power method
    static std::vector<double> vec1(generate_rand(m)); // initial random vector.
    static std::vector<double> vec2(vec1); // initial random vector.

    // initialize solution data
    // initial guess
    std::vector<double> sol(vec1);
    // new solution
    std::vector<double> new_sol(m, 0.);
    // work array
    std::vector<double> work(new_sol);

    // get norm of initial sol
    double alpha    = Tnrm2(m, &sol[0]);
    double gamma    = 1. / alpha;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (mmpi.instancePE0())
        std::cout << "e1:: ITER 0:: = " << alpha << " shft = " << shft
                  << std::endl;

    // residual
    std::vector<double> res(new_sol);
    // initial eigenvalue estimate (for shifted system)
    double beta
        = LinearAlgebraUtils<MemorySpace::Host>::MPdot(m, &sol[0], &new_sol[0]);

    // compute first extent
    int iter1 = 0;
    // loop
    for (int i = 0; i < maxits; i++)
    {
        iter1++;

        // First compute residual for convergence check
        mat.gemv(one, sol, zero, res);
        // store matvec result appropriately scaled for later reuse
        work.clear();
        work.resize(m, 0.);
        LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
            m, gamma, &res[0], &work[0]);
        // Compute residual: res = beta*S*x - mat*x
        (invS_->gramMat())->gemv(beta, sol, -1., res);
        // compute residual norm
        double resnorm = Tnrm2(m, &res[0]);
        // check for convergence
        if (resnorm < 1.0e-2) break;

        // Solve to apply inverse to new_sol to update solution
        // No need to do matvec with scaled copy of sol.
        // Reuse previously stored matvec from residual calculation
        invS_->GramMatLSSolve(&work[0], &new_sol[0]);

        // compute 'shifted' eigenvalue
        beta = LinearAlgebraUtils<MemorySpace::Host>::MPdot(
            m, &sol[0], &new_sol[0]);
        // scale beta by gamma to account for normalizing sol
        beta *= gamma;
        // update solution data
        sol = new_sol;
        // compute norm and update gamma
        alpha = Tnrm2(m, &sol[0]);
        gamma = 1. / alpha;
    }
    // compute first extent (eigenvalue)
    double e1 = beta - shft;
    vec1      = sol;

    // shift matrix by beta and compute second extent
    // store shift
    double shft_e1 = -beta;
    mat.axpy(shft_e1, *(invS_->gramMat()));

    // reset data and begin loop
    sol = vec2;
    new_sol.clear();
    new_sol.resize(m, 0.);
    work  = new_sol;
    alpha = Tnrm2(m, &sol[0]);
    gamma = 1. / alpha;
    beta
        = LinearAlgebraUtils<MemorySpace::Host>::MPdot(m, &sol[0], &new_sol[0]);

    // loop
    if (mmpi.instancePE0())
        std::cout << "e2:: ITER 0:: = " << beta << std::endl;
    int iter2 = 0;
    for (int i = 0; i < maxits; i++)
    {
        iter2++;

        // First compute residual for convergence check
        mat.gemv(one, sol, zero, res);
        // store matvec result appropriately scaled for later reuse
        work.clear();
        work.resize(m, 0.);
        LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
            m, gamma, &res[0], &work[0]);
        // Compute residual: res = beta*S*x - mat*x
        (invS_->gramMat())->gemv(beta, sol, -1., res);
        // compute residual norm
        double resnorm = Tnrm2(m, &res[0]);
        // check for convergence
        if (resnorm < 1.0e-2) break;

        // Solve to apply inverse to new_sol to update solution
        // No need to do matvec with scaled copy of sol.
        // Reuse previously stored matvec from residual calculation
        invS_->GramMatLSSolve(&work[0], &new_sol[0]);

        // compute 'shifted' eigenvalue
        beta = LinearAlgebraUtils<MemorySpace::Host>::MPdot(
            m, &sol[0], &new_sol[0]);
        // scale beta by gamma to account for normalizing sol
        beta *= gamma;
        // update solution data
        sol = new_sol;
        // compute norm and update gamma
        alpha = Tnrm2(m, &sol[0]);
        gamma = 1. / alpha;
    }

    // compute second extent
    double e2 = beta - shft_e1 - shft;
    vec2      = sol;

    // Communicate to take global min and max and save results
    double tmp = e1;
    e1         = std::min(tmp, e2);
    e2         = std::max(tmp, e2);

    mmpi.allreduce(&e1, 1, MPI_MIN);
    mmpi.allreduce(&e2, 1, MPI_MAX);
    //   double tmp = e1;
    //   e1 = min(tmp, e2);
    //   e2 = max(tmp, e2);
    double padding = pad * (e2 - e1);

    if (mmpi.instancePE0())
        std::cout << "Power method Eigen intervals********************  = ( "
                  << e1 << ", " << e2 << ")"
                  << "iter1 = " << iter1 << ", iter2 = " << iter2 << std::endl;

    e1 -= padding;
    e2 += padding;
    interval.push_back(e1);
    interval.push_back(e2);

    // update shft
    shft = std::max(fabs(e1), fabs(e2));

    eig_interval_tm_.stop();
}
