// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SP2.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "ReplicatedMatrix.h"
#include "linear_algebra/blas3_c.h"

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrix.h"
#endif

Timer SP2::getdm_tm_("SP2::getDM");

#ifdef HAVE_BML
double computePartialTrace(bml_matrix_t* A, const std::vector<int>& ids)
{
    assert(A != 0);

    double trace = 0.;
    const int n  = (int)ids.size();
#pragma omp parallel for reduction(+ : trace)
    for (int i = 0; i < n; i++)
    {
        assert(ids[i] < m);
        trace += *((double*)bml_get(A, ids[i], ids[i]));
    }

    return trace;
}
#endif

SP2::SP2(const double tol, const bool distributed)
    : tol_(tol), distributed_(distributed)
{
    // std::cout<<"SP2 with tol = "<<tol_<<endl;
    // Get user defined ratio and tolerance
    Xi_       = nullptr;
    Xi_sq_    = nullptr;
    trace_[0] = 0;
    trace_[1] = 0;
}

SP2::~SP2()
{
    assert(Xi_ != nullptr);
#ifdef HAVE_BML
    bml_deallocate(&Xi_);
    bml_deallocate(&Xi_sq_);
#else
    delete Xi_;
    delete Xi_sq_;
#endif
}

// Calculate A for current Xi_, Xi_sq_
// based on formula in A.M.N. Niklasson, Chapter 16 of
//"Linear scaling techniques in comput. chem. and phys."
// R. Zalesny et al. (eds), 2011
int SP2::calcA(const int nel)
{
    double lhs = std::abs(trace_[1] - (double)nel * 0.5);
    double rhs = std::abs(2. * trace_[0] - trace_[1] - (double)nel * 0.5);
    if (lhs < rhs)
        return 1;
    else
        return 0;
}

void SP2::reduceSumTrace()
{
    if (distributed_)
    {
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.allreduce(&trace_[0], 2, MPI_SUM);
    }
}
// Update Xi_ and Xi_sq_
void SP2::iterate(int A)
{
    if (A)
    {
#ifdef HAVE_BML
        bml_copy(Xi_sq_, Xi_);
#else
        Xi_->copy(*Xi_sq_);
#endif
    }
    else
    {
#ifdef HAVE_BML
        bml_add(Xi_, Xi_sq_, 2., -1., 0.);
#else
        Xi_->scal(2.);
        Xi_->axpy(-1., *Xi_sq_);
#endif
    }

    // update square of density matrix
#ifdef HAVE_BML
    bml_multiply_x2(Xi_, Xi_sq_, 0.);
#else
    Xi_sq_->gemm('n', 'n', 1., *Xi_, *Xi_, 0.);
#endif

    // update traces
#ifdef HAVE_BML
    trace_[0] = computePartialTrace(Xi_, loc_ids_);

    trace_[1] = computePartialTrace(Xi_sq_, loc_ids_);
#else
    trace_[0] = Xi_->computePartialTrace(loc_ids_);

    trace_[1] = Xi_sq_->computePartialTrace(loc_ids_);
#endif
    reduceSumTrace();
}

void SP2::solve(const int nel, const bool verbose)
{
#ifdef PRINT_OPERATIONS
    std::cout << "SP2::solve()..." << std::endl;
#endif

    int i     = 0;
    bool flag = true;

    // main loop
    while (flag && i < 50)
    {
        // Calculate if condition "A" is satisfied
        int A = calcA(nel);

        // Update Xi_ and Xi_sq_
        iterate(A);

        std::cout << std::setprecision(10);
        if (onpe0 && verbose)
            std::cout << "Trace at step " << i << ":" << trace_[0] << std::endl;

        if (fabs(trace_[0] - trace_[1]) < tol_) flag = false;
        i++;
    }
    if (onpe0 && verbose)
    {
        std::cout << "SP2 computed Trace = " << trace_[0] << std::endl;
    }
}

template <>
void SP2::initializeLocalMat(
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& submatM,
    const double emin, const double emax, const std::vector<int>& loc_ids)
{
    loc_ids_ = loc_ids;

    const int n = submatM.n();

#ifdef HAVE_BML
    Xi_    = bml_zero_matrix(dense, double_real, n, n, sequential);
    Xi_sq_ = bml_zero_matrix(dense, double_real, n, n, sequential);
#else
    Xi_       = new SquareLocalMatrices<MATDTYPE, MemorySpace::Host>(1, n);
    Xi_sq_    = new SquareLocalMatrices<MATDTYPE, MemorySpace::Host>(1, n);
#endif

    double factor = 1. / (emax - emin);

    // initialize Xi
    Xi_->copy(submatM);

#ifdef HAVE_BML
    const int n1 = n + 1;
    // shift and scale Xi
    // We shift by -emin since we know it, instead of by the maximum
    // eigenvalue in magnitude as in original paper by Niklasson
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        const int pos = i * n1;
        double* val   = (double*)bml_get(Xi_, i, i);
        *val -= emax;
    }
#else
    Xi_->shift(-emin);
#endif

    // scale
#ifdef HAVE_BML
    bml_scale_inplace(&factor, Xi_);
#else
    Xi_->scal(factor);
#endif

    // initialize Xi_sq_
#ifdef HAVE_BML
    bml_multiply_x2(Xi_, Xi_sq_, 0.);
#else
    (*Xi_sq_).gemm('n', 'n', 1., *Xi_, *Xi_, 0.);
#endif

#ifdef HAVE_BML
    trace_[0] = computePartialTrace(Xi_, loc_ids_);

    trace_[1] = computePartialTrace(Xi_sq_, loc_ids_);
#else
    trace_[0] = Xi_->computePartialTrace(loc_ids_);

    trace_[1] = Xi_sq_->computePartialTrace(loc_ids_);
#endif

    reduceSumTrace();
}

#ifdef MGMOL_USE_SCALAPACK
template <>
void SP2::getDM(dist_matrix::DistMatrix<DISTMATDTYPE>& submatM, // output
    const dist_matrix::DistMatrix<DISTMATDTYPE>& invS)
{
    getdm_tm_.start();

    const int n = invS.n();
    dist_matrix::DistMatrix<DISTMATDTYPE> Xi("Xi", n, n);

    // Here Xi_ is assumed to be "replicated" by each MPI task
    Xi.init(Xi_->getSubMatrix(), n);

    submatM.gemm('n', 'n', 2., Xi, invS, 0.);

    getdm_tm_.stop();
}
#endif

template <>
void SP2::getDM(ReplicatedMatrix& submatM, // output
    const ReplicatedMatrix& invS)
{
    getdm_tm_.start();

    const int n = Xi_->n();
    ReplicatedMatrix Xi("Xi", n, n);

    Xi.init(Xi_->getSubMatrix(), n);

    submatM.gemm('n', 'n', 2., Xi, invS, 0.);

    getdm_tm_.stop();
}
