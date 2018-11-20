// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ProjectedMatricesSparse.h"
#include <cassert>
#include <fstream>
#include <iomanip>

#include "BasicDataDistributors.h"
#include "Control.h"
#include "HDFrestart.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
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

ProjectedMatricesSparse::ProjectedMatricesSparse(
    const int ndim, LocalizationRegions* lrs, ClusterOrbitals* local_cluster)
{
    assert(lrs != 0);

    dim_     = ndim;
    min_val_ = 0.25;

    lrs_              = lrs;
    local_cluster_    = local_cluster;
    distributor_matS_ = BasicDataDistributors::gramMatrixDistributor();
    distributor_invS_
        = BasicDataDistributors::centeredOrbitalsOverlapDistributor();
    distributor_sH_ = BasicDataDistributors::orbitalsProdWithHDistributor();

    invS_ = 0;

    dm_ = 0;

    sH_      = 0;
    matHB_   = 0;
    submatT_ = 0;
    localX_  = 0;
    localT_  = 0;

    isDataSetup_ = false;
}

void ProjectedMatricesSparse::clearData()
{
    //    if(isDataSetup_)
    //    {
    // if(onpe0)cout<<"delete invS"<<endl;
    //       delete invS_; invS_=0;
    //       delete dm_; dm_ = 0;
    // if(onpe0)cout<<"delete localX"<<endl;
    delete localX_;
    localX_ = 0;
    // if(onpe0)cout<<"delete localT"<<endl;
    delete localT_;
    localT_ = 0;
    // if(onpe0)cout<<"delete sH"<<endl;
    delete sH_;
    sH_ = 0;
    // if(onpe0)cout<<"delete HB"<<endl;
    delete matHB_;
    matHB_ = 0;
    // if(onpe0)cout<<"delete T"<<endl;
    delete submatT_;
    submatT_ = 0;

    isDataSetup_ = false;
    //    }
}

ProjectedMatricesSparse::~ProjectedMatricesSparse()
{
    assert(invS_ != 0);
    assert(dm_ != 0);
    assert(localX_ != 0);
    assert(localT_ != 0);
    assert(submatT_ != 0);
    assert(sH_ != 0);
    assert(matHB_ != 0);
    assert(distributor_sH_ != 0);
    assert(distributor_matS_ != 0);
    assert(distributor_invS_ != 0);

    // if(onpe0)cout<<"delete invS"<<endl;
    delete invS_;
    invS_ = 0;
    delete dm_;
    dm_ = 0;
    clearData();
};

void ProjectedMatricesSparse::setup(
    const double kbt, const int nel, const vector<vector<int>>& global_indexes)
{
    // assert( (short)global_indexes.size()>0 );
    // assert( (short)global_indexes[0].size()>0 );
    Control& ct = *(Control::instance());

    ProjectedMatricesInterface::setup(kbt, nel, global_indexes);

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
        invS_ = 0;
        delete dm_;
        dm_ = 0;
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

        lsize_ = locvars_.size(); // table.get_size();

        // reset invS and DM data
        if (invS_ == 0 || dm_ == 0 || ct.AtomsDynamic() == AtomsDynamicType::MD)
        {
            invS_ = new ShortSightedInverse((*lrs_), locvars_, local_cluster_);
            dm_   = new DensityMatrixSparse(
                (*lrs_), dim_, locvars_, local_cluster_);
        }

        matHB_ = new VariableSizeMatrix<sparserow>("HB", lsize_);
        // estimate size of table needed for efficient access to elements of sH
        sH_ = new VariableSizeMatrix<sparserow>("sH", 4096);
        /* initialize Sparse H matrix -- this is necessary for efficient data
         * distribution */
        (*sH_).setupSparseRows(locvars_);
        //         (*sH_).setupsparserows(locfcns);

        submatT_ = new VariableSizeMatrix<sparserow>("Theta", lsize_);

        localX_ = new SquareLocalMatrices<MATDTYPE>(subdiv_, chromatic_number_);
        localT_ = new SquareLocalMatrices<MATDTYPE>(subdiv_, chromatic_number_);
    }
    // get min and max of functions centered on local processors
    ///*
    int locmin, locmax;
    locmin = locfcns.size();
    locmax = locmin;
    mmpi.allreduce(&locmin, 1, MPI_MIN);
    mmpi.allreduce(&locmax, 1, MPI_MAX);

    if (onpe0 && ct.verbose > 0)
    {
        printf("Max. number of locally centered functions: %d \n", locmax);
        printf("Min. number of locally centered functions: %d \n", locmin);
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
    assert(invS_ != NULL);

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

    return;
}
double ProjectedMatricesSparse::getExpectationH()
{
    assert(invS_ != NULL);
    assert(matHB_ != 0);
    //   assert(matHB_->n() != 0);

    return dm_->getTraceDotProductWithMat(matHB_);
}

double ProjectedMatricesSparse::getExpectationMat(
    VariableSizeMatrix<sparserow>* mat)
{
    assert(invS_ != NULL);
    assert(mat != 0);
    assert(mat->n() != 0);

    return dm_->getTraceDotProductWithMat(mat);
}

void ProjectedMatricesSparse::consolidateOrbitalsOverlapMat(
    VariableSizeMatrix<sparserow>& mat)
{
    consolidate_H_tm_.start();
    std::vector<int> locfcns;
    (*lrs_).getLocalSubdomainIndices(locfcns);

    // update locally centered row data
    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    // consolidate locally centered data
    DataDistribution distributor1(
        "overlap", (*lrs_).max_radii(), myPEenv, domain);
    distributor1.augmentLocalData(mat, false);

    // now gather updated data from neighbors.
    // gather from neighbors that contain functions that overlap with
    // locally centered functions
    Control& ct(*Control::instance());
    DataDistribution distributor("overlap2", ct.spread_radius, myPEenv, domain);

    mat.consolidate(locfcns, distributor);

    consolidate_H_tm_.stop();
}
void ProjectedMatricesSparse::consolidateH()
{
    /* Gather data to initialize matH and matHB */

    assert(sH_->n() == (int)locvars_.size());

    consolidate_H_tm_.start();
    std::vector<int> locfcns;
    (*lrs_).getLocalSubdomainIndices(locfcns);

    /* gather data for matH amd matHB */
    (*distributor_sH_).augmentLocalData((*sH_), false);

    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    // now gather updated data from neighbors.
    // gather from neighbors that contain functions that overlap with locally
    // centered functions
    Control& ct(*(Control::instance()));
    DataDistribution distributor("H", ct.spread_radius, myPEenv, domain);

    sH_->consolidate(locfcns, distributor);

    (*matHB_).copyData((*sH_), lsize_);

    consolidate_H_tm_.stop();
    return;
}

// computes the trace of the dotproduct with invS
// assumes matrix ss contains partial contributions
double ProjectedMatricesSparse::dotProductWithInvS(
    const SquareLocalMatrices<MATDTYPE>& local_product)
{
    return computeTraceInvSmultMat(local_product, true);
}

// computes the trace of the dotproduct with DM
// assumes matrix ss contains partial contributions
double ProjectedMatricesSparse::dotProductWithDM(
    const SquareLocalMatrices<MATDTYPE>& local_product)
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
    const SquareLocalMatrices<MATDTYPE>& local_product)
{
    std::cerr << "ERROR: ProjectedMatricesSparse::dotProductSimple() \
               not implemented!!!"
              << std::endl;

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Abort(mmpi.commSameSpin(), 0);

    return -1.;
}

void ProjectedMatricesSparse::printMatrices(ostream& os) const
{
    printS(os);
    printH(os);
    printTheta(os);
    printHB(os);

    /* print stats for Gram Matrix data distribution */
    Control& ct = *(Control::instance());
    if ((ct.verbose > 1) && onpe0)
    {
        cout << " Gram Matrix data distribution stats " << endl;
        (*distributor_matS_).printStats();
    }
}

double ProjectedMatricesSparse::getNel() const { return 2. * dim_; }

double ProjectedMatricesSparse::getEigSum()
{
    if (dim_ == 0) return 0.;
    assert(invS_ != NULL);
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
    assert((*invS_).gramMat() != 0);

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
        cout << "" << endl;
        cout << "ProjectedMatricesSparse::getLinDependent2states(): MPI task "
             << myid << endl;
        cout << "eigenvector size: " << evec.size() << endl;
        // cout<<"locst1="<<locst1<<", locst2="<<locst2<<endl;
        cout << "off-diagonal element: " << vmax2 << endl;
        cout << "2x2 submatrix associated with small eigenvalue:" << endl;
        mat->printMatBlock2(st1, st2, cout);
        vector<double> tmp(bsize);
        LSMat.matvec(&evec[0], &tmp[0]);
        int one      = 1;
        double norme = dnrm2(&bsize, &evec[0], &one);
        double norma = dnrm2(&bsize, &tmp[0], &one);
        cout << "eigenvalue: " << out.val << endl;
        cout << "Norm S*v: " << norma / norme << endl;
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

void ProjectedMatricesSparse::printTimers(ostream& os)
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
    SquareLocalMatrices<MATDTYPE>* localM)
{
    // return if matrix size is zero
    if (localM->n() == 0) return;

    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        MATDTYPE* localM_iloc = localM->getSubMatrix(iloc);
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
    const SquareLocalMatrices<MATDTYPE>& mat, const bool consolidate)
{
    Control& ct = *(Control::instance());
    VariableSizeMatrix<sparserow> vsmat("VS", lsize_);
    vsmat.setupSparseRows(locvars_);
    vsmat.initializeMatrixElements(mat, global_indexes_, ct.numst);

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
    const SquareLocalMatrices<MATDTYPE>& mat, const bool consolidate)
{
    Control& ct = *(Control::instance());
    VariableSizeMatrix<sparserow> vsmat("VS", lsize_);
    vsmat.setupSparseRows(locvars_);
    vsmat.initializeMatrixElements(mat, global_indexes_, ct.numst);

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

double ProjectedMatricesSparse::computeTraceInvSmultMatmultTheta(
    SquareLocalMatrices<MATDTYPE>& mat, const bool consolidate)
{
    /*
       assert(localT_ != 0);
       assert(submatT_ != 0);

       double trace = 0.;

       if(localT_ != 0)
       {
          // compute prod = mat * theta
          SquareLocalMatrices<MATDTYPE> prod(subdiv_, chromatic_number_);
          prod.gemm('n', 'n', 1.0, mat, *localT_, 0.);

          // compute invS * prod
          VariableSizeMatrix<sparserow> vsmat(lsize_);
          vsmat.setupSparseRows(locvars_);
          vsmat.initializeMatrixElements(prod, global_indexes_);

          if(consolidate)
          {
             //consolidate only locally centered data since trace is needed
             DataDistribution distributor((*lrs_).max_radii());
             distributor.augmentLocalData(vsmat, false);
          }
       }



       std::vector<int>locfcns;
       (*lrs_).getLocalSubdomainIndices(locfcns);
       // try initializing with locfcns
       VariableSizeMatrix<sparserow> Cmat(lsize_);
       Cmat.setupSparseRows(locvars_);
       // do product of mat with theta
       // get pointer to new Gram matrix
    //   VariableSizeMatrix<sparserowtab> *gmatptr = (*invS_).gramMat();
    //   vsmat.AmultSymBLocal(submatT_, Cmat, locfcns, *gmatptr, false);

       // compute trace with invS
    //   double trace = 0.;


       return getExpectationMat(&Cmat);
    */
    return 0.;
}

void ProjectedMatricesSparse::applyInvS(SquareLocalMatrices<MATDTYPE>& mat)
{
    Control& ct = *(Control::instance());
    mat.transpose();

    // convert mat into variablesizematrix object
    VariableSizeMatrix<sparserow> vsmat("vs", lsize_);
    vsmat.setupSparseRows(locvars_);
    vsmat.initializeMatrixElements(mat, global_indexes_, ct.numst);

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

    return;
}

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
    assert(invS_ != 0);
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
    double alpha = Tnrm2(m, &sol[0]);
    double gamma = 1. / alpha;
    if (onpe0)
        cout << "e1:: ITER 0:: = " << alpha << " shft = " << shft << endl;

    // residual
    std::vector<double> res(new_sol);
    // initial eigenvalue estimate (for shifted system)
    double beta = MPdot(m, &sol[0], &new_sol[0]);

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
        MPaxpy(m, gamma, &res[0], &work[0]);
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
        beta = MPdot(m, &sol[0], &new_sol[0]);
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
    beta  = MPdot(m, &sol[0], &new_sol[0]);

    // loop
    if (onpe0) cout << "e2:: ITER 0:: = " << beta << endl;
    int iter2 = 0;
    for (int i = 0; i < maxits; i++)
    {
        iter2++;

        // First compute residual for convergence check
        mat.gemv(one, sol, zero, res);
        // store matvec result appropriately scaled for later reuse
        work.clear();
        work.resize(m, 0.);
        MPaxpy(m, gamma, &res[0], &work[0]);
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
        beta = MPdot(m, &sol[0], &new_sol[0]);
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
    e1         = min(tmp, e2);
    e2         = max(tmp, e2);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&e1, 1, MPI_MIN);
    mmpi.allreduce(&e2, 1, MPI_MAX);
    //   double tmp = e1;
    //   e1 = min(tmp, e2);
    //   e2 = max(tmp, e2);
    double padding = pad * (e2 - e1);

    if (onpe0)
        cout << "Power method Eigen intervals********************  = ( " << e1
             << ", " << e2 << ")"
             << "iter1 = " << iter1 << ", iter2 = " << iter2 << endl;

    e1 -= padding;
    e2 += padding;
    interval.push_back(e1);
    interval.push_back(e2);

    // update shft
    shft = max(fabs(e1), fabs(e2));

    eig_interval_tm_.stop();
}
