// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "Control.h"
#include "LinearSolverMatrix.h"
#include "MGmol_MPI.h"
#include "MGmol_blas1.h"
#include "PreconILU.h"
#include "ShortSightedInverse.h"
#include "Table.h"
#include "VariableSizeMatrix.h"
#include "tools.h"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

#define epsmac 1.0e-16

Timer ShortSightedInverse::Gram_Matrix_data_distribution_tm_(
    "ShortSightedInverse::GramMat_data_distribution");
Timer ShortSightedInverse::reset_tm_("ShortSightedInverse::reset");
Timer ShortSightedInverse::inverse_solve_tm_("ShortSightedInverse::solve");
Timer ShortSightedInverse::compute_invS_tm_(
    "ShortSightedInverse::Compute_Inverse");
Timer ShortSightedInverse::gather_invS_tm_(
    "ShortSightedInverse::gather_inverse");
Timer ShortSightedInverse::gmres_tm_("ShortSightedInverse::gmres");
Timer ShortSightedInverse::pc_setup_tm_("ShortSightedInverse::pc_setup");
// Timer   ShortSightedInverse::pc_solve_tm_("ShortSightedInverse::pc_solve");
Timer ShortSightedInverse::linear_solver_matvec_tm_(
    "ShortSightedInverse::linear_solver_matvec");
Timer ShortSightedInverse::linear_solver_matrix_init_tm_(
    "ShortSightedInverse::linear_solver_matrix_init");

// const double mat_tol = 1.0e-14;
ShortSightedInverse::ShortSightedInverse(
    std::shared_ptr<LocalizationRegions> lrs, const std::vector<int>& locvars,
    ClusterOrbitals* local_cluster)
    : locvars_(locvars)
{
    loc_radius_ = lrs->max_radii();
    max_size_   = MAX_MAT_SIZE;

    /* preconditioner and solver parameters */
    Control& ct = *(Control::instance());

    droptol_    = (float)ct.ilu_droptol;
    MaxFil_     = ct.ilu_maxfil;
    lof_        = ct.ilu_lof;
    fgmres_tol_ = ct.fgmres_tol;
    im_         = ct.fgmres_kim;
    maxits_     = ct.fgmres_maxits;
    if (ct.ilu_type == 0)
        ilutype_ = PCILUK;
    else if (ct.ilu_type == 1)
        ilutype_ = PCMILUT;
    else if (ct.ilu_type == 2)
        ilutype_ = PCILUT;
    else
        ilutype_ = PCILU0;

    /* set local data size */
    lsize_ = (int)locvars.size();
    /* Get local functions centered on this processor */
    if (local_cluster != nullptr)
        locfcns_ = local_cluster->getClusterIndices();
    else
        lrs->getLocalSubdomainIndices(locfcns_);

    /* Setup/ Initialize some local objects */
    /* setup Gram matrix */
    if (onpe0 && ct.verbose > 0) std::cout << "lsize_=" << lsize_ << std::endl;
    // estimate of size parameter (4096) needed for efficient table element
    // access
    gramMat_ = new VariableSizeMatrix<sparserowtab>("Gram", 4096);
    /* setup invS matrix for storing inverse */
    invS_  = new VariableSizeMatrix<sparserow>("invS", lsize_);
    matLS_ = new LinearSolverMatrix<lsdatatype>(0, 0);

    issetup_        = true;
    isInvSUpToDate_ = false;

    resnorm_ = 0.0;
    precon_  = nullptr;
    recompute_pc_
        = true; /* recompute the preconditioner or not - initially true */
}

ShortSightedInverse::~ShortSightedInverse()
{
    delete gramMat_;
    gramMat_ = nullptr;
    delete invS_;
    invS_ = nullptr;
    delete matLS_;
    matLS_ = nullptr;
    delete precon_;
    precon_ = nullptr;
}

/* Reset/Setup local objects */
void ShortSightedInverse::reset()
{
    reset_tm_.start();
    /* reset Gram matrix */
    (*gramMat_).reset();

    /* reset preconditioner */
    if (recompute_pc_ == true)
    {
        if (precon_ != nullptr)
        {
            delete precon_;
            precon_ = nullptr;
        }
    }

    reset_tm_.stop();
    return;
}

// int ShortSightedInverse::solve(LinearSolverMatrix<lsdatatype>& LSMat)
int ShortSightedInverse::solve()
{
    inverse_solve_tm_.start();

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Control& ct     = *(Control::instance());

    int conv       = 0;
    int maxits     = 0;
    double maxrnrm = 0.;
    std::vector<int> its(locfcns_.size());
    std::vector<double> rnrm(locfcns_.size());

    const unsigned locfcns_size = locfcns_.size();
#pragma omp parallel
    {
        // create Linear solver object
        LinearSolver solver;

        const int m    = aug_size_;
        double* solptr = new double[m];

#pragma omp for reduction(+ : conv)
        for (unsigned int i = 0; i < locfcns_size; i++)
        {
            int* rindex       = (int*)(*gramMat_).getTableValue(locfcns_[i]);
            const int lrindex = *rindex;

            /* call gmres */
            gmres_tm_.start();

            memset(solptr, 0, m * sizeof(double));
            conv += solver.solve((*matLS_), (*precon_), lrindex, solptr,
                fgmres_tol_, im_, maxits_);

            its[i] = solver.iters();

            rnrm[i] = solver.residualNorm();
            gmres_tm_.stop();

            /* Update invS */
            const int* row = (int*)(*invS_).getTableValue(locfcns_[i]);
            assert(row != nullptr);
            const int* cols = (*gramMat_).rowIndexes();
            (*invS_).initializeLocalRow(m, *row, cols, solptr);
        }
        delete[] solptr;

    } // end OpenMP region

    for (unsigned int i = 0; i < locfcns_size; i++)
    {
        maxits = (maxits > its[i]) ? maxits : its[i];
    }
    for (unsigned int i = 0; i < locfcns_size; i++)
    {
        maxrnrm = (maxrnrm > rnrm[i]) ? maxrnrm : rnrm[i];
    }

    /* check for convergence */

    int tmp[2] = { maxits, conv };
    mmpi.allreduce(&tmp[0], 2, MPI_MAX);
    maxits = tmp[0];
    conv   = tmp[1];

    mmpi.allreduce(&maxrnrm, 1, MPI_MAX);

    inverse_solve_tm_.stop();

    // check for convergence
    if (conv == 0)
    {
        recompute_pc_ = false;
        if ((onpe0) && (ct.verbose > 1))
            printf("GMRES converged in %d iterations. Max .Final res. norm = "
                   "%2.8e \n",
                maxits, maxrnrm);
        return 0;
    }

    // unconverged solution
    if (recompute_pc_ == false)
    {
        delete precon_;
        precon_       = nullptr;
        recompute_pc_ = true;
        if (onpe0)
            printf(
                "WARNING: Recompute Preconditioner - GMRES did not converge "
                "after max. its =  %d iterations. Final res. norm = %2.8e \n",
                maxits_, maxrnrm);
        return 1;
    }
    else
    {
        if (onpe0)
            printf("WARNING: GMRES did not converge after max. its =  %d "
                   "iterations. Final res. norm = %2.8e \n",
                maxits_, maxrnrm);
        return 0;
    }
}

/* compute the inverse of gramMat_ */
void ShortSightedInverse::computeInvS(DataDistribution& distributor_invS)
{
    assert(issetup_ == true);
    if (isGramAugmented_ == true) assert((*gramMat_).nzmax() > 0);

    if (isInvSUpToDate_) return;

    compute_invS_tm_.start();

    Control& ct = *(Control::instance());

    new_pc_ = false;

    /* perform data distribution */
    if (isGramAugmented_ == false)
    {
        augmentGramMatrix();
    }

    /* prepare invS_ for inverse data - match row pattern of local gram matrix
     */

    /* setup sparse local inverse matrix */
    (*invS_).setupSparseRows(locvars_);

    /* setup linear solver matrix */
    matLS_->reset(aug_size_, (*gramMat_).nnzmat());
    //    LinearSolverMatrix<lsdatatype> LSMat(aug_size_, (*gramMat_).nnzmat());

    /* initialize linear solver matrix */
    linear_solver_matrix_init_tm_.start();
    //    LSMat.init((*gramMat_), false);
    matLS_->init((*gramMat_), false);
    linear_solver_matrix_init_tm_.stop();

    /*** dump the matrix for debug **/
    /*
    if(onpe0){
    FILE *ft = fopen("mat0.mm", "w");
    fprintf(ft," %d %d %d \n", LSMat.n(), LSMat.n(), LSMat.nnz());
    for(int i=0; i<LSMat.n(); i++)
    {
       bool flag = false;
       for(int k=LSMat.nzptrval(i); k<LSMat.nzptrval(i+1); k++)
       {
          if(i == LSMat.getColumnIndex(k))
              flag = true;
          fprintf(ft," %d %d %2.8e \n", i, LSMat.getColumnIndex(k),
    LSMat.getRowEntry(k));
       }
       if(flag == false)
          fprintf(ft," %d %d %2.8e \n", i, i, MAT_TOL);
    }
    fclose(ft);
    }
    MPI_Finalize();
    exit(0);
    */
    /* update preconditioner */

    /* update inverse - solve linear system */

    int solve_flag = 1;
    while (solve_flag)
    {
        if (recompute_pc_ == true)
        {
            pc_setup_tm_.start();
            assert(precon_ == nullptr);
            //          precon_ = new
            //          PreconILU<pcdatatype>(LSMat,droptol_,MaxFil_,lof_);
            precon_
                = new PreconILU<pcdatatype>((*matLS_), droptol_, MaxFil_, lof_);
            /* setup the preconditioner */
            //          (*precon_).setup(LSMat, ilutype_);
            (*precon_).setup((*matLS_), ilutype_);
            if (onpe0 && (ct.verbose > 1))
                (*MPIdata::sout)
                    << std::scientific << "ILU preconditioner fill_factor = "
                    << ((double)(*precon_).nnz_ilu()) / ((double)matLS_->nnz())
                    << std::endl;
            new_pc_ = true;
            pc_setup_tm_.stop();
        }

        //        solve_flag = solve((*matLS_));
        solve_flag = solve();

    } // end while

    /* gather inverse data */
    gather(distributor_invS);
    isInvSUpToDate_ = true;

    compute_invS_tm_.stop();

    return;
}

/* initialize the gram matrix */
void ShortSightedInverse::initGramMatrix(
    const LocalMatrices<MATDTYPE, MemorySpace::Host>& ss,
    const std::vector<std::vector<int>>& global_indexes,
    const int new_orbitals_index)
{
    if ((*gramMat_).n() > 0 && gm_orbitals_index_ == new_orbitals_index) return;

    /* reset the local data */
    reset();

    /* update Gram Matrix */
    Control& ct = *(Control::instance());
    (*gramMat_).setupSparseRows(locvars_);
    (*gramMat_).insertMatrixElements(
        ss, global_indexes, ct.numst); // should lsize_ be an argument here??

    /* set augmented gram matrix to false */
    isGramAugmented_ = false;

    /* set inverse to be out-dated */
    isInvSUpToDate_ = false;

    /* update orbitals_index_ */
    gm_orbitals_index_ = new_orbitals_index;
}

/* gather inverse contributions for local Gram matrix */
void ShortSightedInverse::gather(DataDistribution& distributor_invS)
{
    assert(issetup_ == true);
    gather_invS_tm_.start();
    /* Distribute local inverse and gather data */
    distributor_invS.updateLocalRows((*invS_));
    gather_invS_tm_.stop();
}

/* perform the trace of the product of some diagonal matrix represented by the
 * vector ddiag with the diagonal of invS. It is assumed that ddiag contains the
 * diagonal entries of the functions centered on the local processor.(ie.
 * functions in locfcns_)
 */
double ShortSightedInverse::getTraceDiagProductWithInvS(
    std::vector<DMDTYPE>& ddiag)
{
    assert((int)ddiag.size() == (int)locfcns_.size());

    double trace = 0.0;
    std::vector<int>::iterator itr;
    int pos;
    for (itr = locfcns_.begin(), pos = 0; itr != locfcns_.end(); ++itr, pos++)
    {
        trace += ddiag[pos] * (*invS_).get_value(*itr, *itr);
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();

    double val = 0.0;
    MPI_Allreduce(&trace, &val, 1, MPI_DOUBLE, MPI_SUM, comm);

    return val;
}

/* print a few rows of the Gram matrix */
void ShortSightedInverse::printGramMat(std::ostream& os, int nrows) const
{
    if (isGramAugmented_ == true)
    {
        if (onpe0)
            os << " Local GramMatrix ... n = " << lsize_
               << " augmented_size = " << aug_size_ << std::endl;
        (*gramMat_).print(os, locfcns_, nrows);
    }
    else
    {
        const int n = (*gramMat_).n();
        /* augment gram matrix copy and print */
        /* only need to gather data for functions centered on local PE */
        VariableSizeMatrix<sparserow> lmat((*gramMat_), true);
        /* define data distribution object */
        Mesh* mymesh             = Mesh::instance();
        const pb::Grid& mygrid   = mymesh->grid();
        const pb::PEenv& myPEenv = mymesh->peenv();
        double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
        DataDistribution distributor("printGram", loc_radius_, myPEenv, domain);

        /* gather data -- need to update only centered functions so append=false
         */
        bool append
            = false; /* gather in one dimension only - assumes symmetry */
        distributor.augmentLocalData(lmat, append);

        if (onpe0) os << " Local GramMatrix size ... n = " << n << std::endl;

        lmat.print(os, locfcns_, nrows);
    }
    return;
}

double ShortSightedInverse::getTraceDotProductWithInvS(
    VariableSizeMatrix<sparserow>* mat)
{
    /* compute trace */
    double trace = 0.0;
    for (std::vector<int>::iterator itr = locfcns_.begin();
         itr != locfcns_.end(); ++itr)
    {
        trace += (*invS_).AmultSymB_ij(mat, *itr, *itr);
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();

    double val = 0.0;
    MPI_Allreduce(&trace, &val, 1, MPI_DOUBLE, MPI_SUM, comm);

    return val;
}

void ShortSightedInverse::augmentGramMatrix()
{
    Gram_Matrix_data_distribution_tm_.start();

    // Alternate approach:
    // Gather to complete only locally centered rows/cols
    // Distribute/ receive complete rows from neighbors within spread_radius.
    // This latter stage should only do row/col assignments (no inserts of
    // individual entries)
    // TODO: Consider setting up sparse rows for needed data prior to second
    // stage communication.
    //       This should speed it up a bit.

    Control& ct = *(Control::instance());
    // update locally centered row data
    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    // consolidate locally centered data
    DataDistribution distributor1("overlap", loc_radius_, myPEenv, domain);
    distributor1.augmentLocalData((*gramMat_), false);

    const bool print_flag = (onpe0 && ct.verbose > 2);
    if (print_flag)
    {
        std::cout << " Gram Matrix data distribution stats: distributor1 "
                  << std::endl;
        distributor1.printStats();
    }

    DataDistribution distributor2(
        "overlap2", ct.spread_radius, myPEenv, domain);
    distributor2.consolidateMatrix(locfcns_, (*gramMat_));

    if (print_flag)
    {
        std::cout << " Gram Matrix data distribution stats: distributor2 "
                  << std::endl;
        distributor2.printStats();
    }

    isGramAugmented_ = true;
    aug_size_        = (*gramMat_).n();

    Gram_Matrix_data_distribution_tm_.stop();
}

void ShortSightedInverse::printTimers(std::ostream& os)
{
    compute_invS_tm_.print(os);
    inverse_solve_tm_.print(os);
    gather_invS_tm_.print(os);
    gmres_tm_.print(os);
    pc_setup_tm_.print(os);
    //   pc_solve_tm_.print(os);
    linear_solver_matvec_tm_.print(os);
    linear_solver_matrix_init_tm_.print(os);
    Gram_Matrix_data_distribution_tm_.print(os);
    reset_tm_.print(os);
}

void ShortSightedInverse::printGramMM(std::ofstream& tfile)
{
    assert(gramMat_ != nullptr);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int n           = locfcns_.size();
    mmpi.allreduce(&n, 1, MPI_SUM);

    /* augment gram matrix copy and print */
    /* only need to gather data for functions centered on local PE */
    VariableSizeMatrix<sparserow> lmat((*gramMat_), true);
    /* define data distribution object */
    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    DataDistribution distributor("GramMM", loc_radius_, myPEenv, domain);

    /* gather data */
    bool append = false; /* gather in one dimension only - assumes symmetry */
    distributor.augmentLocalData(lmat, append);

    int nzmat = 0;
    for (std::vector<int>::const_iterator it = locfcns_.begin();
         it != locfcns_.end(); ++it)
    {
        int* row = (int*)lmat.getTableValue(*it);
        nzmat += lmat.nnzrow(*row);
    }
    mmpi.allreduce(&nzmat, 1, MPI_SUM);

    if (onpe0) tfile << n << " " << n << " " << nzmat << std::endl;

    lmat.print(tfile, locfcns_, n);

    return;
}

/* compute condition number of gramMat_ */
double ShortSightedInverse::computeGramMatCond()
{
    assert(isGramAugmented_ == true);

    double condest = 0.;
    /* get size of submatrix and set its extents */
    const int bsize  = locvars_.size();
    const int bnzmax = (*gramMat_).getNzmaxSubmat(0, bsize - 1);
    /* construct Linear Solver matrix from gram matrix*/

    LinearSolverMatrix<lsdatatype> LSMat(bsize, bnzmax);
    LSMat.initSquareMat((*gramMat_), false);
    // Linear solver object to compute condition number
    LinearSolver solver;
    // compute "local" condition number (estimate)
    if (LSMat.n() > 0) condest = solver.computeConditionNumber(LSMat);

    // get max condition number
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&condest, 1, MPI_MAX);
    return condest;
}

std::vector<int> ShortSightedInverse::centeredFcnLocalIds()
{
    assert(gramMat_->n() != 0);

    std::vector<int> local_cols;
    for (unsigned int i = 0; i < locfcns_.size(); i++)
    {
        int* rindex = (int*)(*gramMat_).getTableValue(locfcns_[i]);
        local_cols.push_back(*rindex);
    }

    return local_cols;
}

// Do local solve with (augmented) Gram Matrix
int ShortSightedInverse::GramMatLSSolve(const double* rhs, double* sol)
{
    Control& ct = *(Control::instance());

    int conv    = 0;
    int its     = 0;
    double rnrm = 0.;

    // create Linear solver object
    LinearSolver solver;

    const int m = matLS_->n();
    /* call gmres */
    gmres_tm_.start();

    memset(sol, 0, m * sizeof(double));
    conv = solver.solve(
        (*matLS_), (*precon_), rhs, sol, fgmres_tol_, im_, maxits_);
    its  = solver.iters();
    rnrm = solver.residualNorm();
    gmres_tm_.stop();

    /* check for convergence */

    if (conv == 0)
    {
        if ((onpe0) && (ct.verbose > 3))
            printf("GMRES converged in %d iterations. Max .Final res. norm = "
                   "%2.8e \n",
                its, rnrm);
        return 0;
    }
    else // unconverged solution
    {
        if ((onpe0) && (ct.verbose > 3))
            printf("WARNING: GMRES did not converge after max. its =  %d "
                   "iterations. Final res. norm = %2.8e \n",
                maxits_, rnrm);
        return 0;
    }
}
