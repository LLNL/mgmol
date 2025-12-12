// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cassert>
#include <iomanip>
#include <iostream>

#include "MGmol_MPI.h"
#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "NOLMOTransform.h"
#include "mputils.h"

// #define DEBUG 1

////////////////////////////////////////////////////////////////////////////////

void NOLMOTransform::init_transform(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& ls, const bool reset_positions)
{
    if (bcr_->onpe0())
        (*MPIdata::sout) << "NOLMOTransform::init_transform()" << std::endl;

    const dist_matrix::DistMatrix<DISTMATDTYPE>& work(ls);
    dist_matrix::DistMatrix<DISTMATDTYPE> tmp(ls);
    tmp.transpose(1., work, 0.);
    tmp.trset('u'); // set to zero lower part of upper triangular matrix

    tmp.allgather(&lst_[0], nst_);

    // copy lst_ into a_
    *a_ = tmp;

    if (bcr_->mype() < a_->npcol())
        if (reset_positions || !set_positions_)
        {
            int n = (int)ara0_[0].size();
            if (bcr_->onpe0()) (*MPIdata::sout) << "set ara0..." << std::endl;
            std::vector<double> alphan0(n);
            std::vector<double> inv_alphan0(n);
            getNorms(alphan0, inv_alphan0);

            for (int k = 0; k < 2 * NDIM; k++)
            {
                for (int i = 0; i < n; i++)
                {
                    const int ii = i * nst_ + i + offset_;
                    ara0_[k][i]  = r_[k]->val(ii) * inv_alphan0[i];
                }
            }
            set_positions_ = true;
        }
}

////////////////////////////////////////////////////////////////////////////////
double NOLMOTransform::spread2(int i, int j) const
{
    assert(i >= 0 && i < nst_);
    assert(j >= 0 && j < NDIM);
    assert(cell_[j] > 0.);

    int ii = -1;

    // integral over domain of position operator^2 (sin^2+cos^2) * 0.5 * M_1_PI
    const double lby2pi = 0.5 * M_1_PI * cell_[j];

    double cs[2] = { 0., 0. };
    double inv_norm;
    if (i >= offset_ && i < offset_ + lnst_ && offset_ >= 0)
    {
        const double alphan = a_->dotColumns(i - offset_, i - offset_);
        assert(alphan > 0.);
        inv_norm = 1. / alphan;

        ii = nst_ * (i - offset_) + i;
        assert(ii >= 0);
        cs[0] = r_[2 * j]->val(ii) * inv_norm;
        cs[1] = r_[2 * j + 1]->val(ii) * inv_norm;
    }

#ifdef DEBUG
    if (ii != -1)
    {
        double c2 = b_[2 * j]->val(ii) * inv_norm;
        double s2 = b_[2 * j + 1]->val(ii) * inv_norm;
        (*MPIdata::sout) << "s2+c2=" << s2 + c2 << std::endl;
        (*MPIdata::sout) << "lby2pi*lby2pi=" << lby2pi * lby2pi << std::endl;
    }
#endif

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (offset_ > 0 && ii != -1)
    {
        mmpi.send(&cs[0], 2, 0);
    }
    if (offset_ == 0 && ii == -1)
    {
        int src = i / bsize_;
        mmpi.recv(&cs[0], 2, src);
    }
    if (bcr_->npcol() > 1) mmpi.bcast(cs, 2);

    return lby2pi * lby2pi - (cs[0] * cs[0] + cs[1] * cs[1]);
    // return s2+c2 - (c*c + s*s);
}

////////////////////////////////////////////////////////////////////////////////
void NOLMOTransform::compute_transform(const int maxsweep, const double tol)
{
    // sine and cosine operators only: compute NOLMO transform

    double delta;
    // delta= nolmo(maxsweep,tol);
    delta = nolmo_fixedCenters(maxsweep, tol);

    if (bcr_->onpe0())
        if (delta > tol)
            (*MPIdata::sout)
                << " NOLMOTransform: decrease was " << std::setprecision(8)
                << delta << " after " << maxsweep << " iterations" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

double NOLMOTransform::nolmo(int maxsweep, double tol)
{
    // The input matrices are given in r_[k] and b_[k], k = 0,..,m-1
    // the initial non-orthogonal transformation is stored in a_
    // the non-orthogonal transformation is stored in mat_[i+nst_*j]

    const int m = 2 * NDIM;

    DISTMATDTYPE one = 1.;

    std::vector<double> alphan(lnst_);
    std::vector<double> inv_alphan(lnst_);

    dist_matrix::DistMatrix<DISTMATDTYPE> grads(
        "grad", *bcr_, nst_, nst_, nst_, bsize_);
    const int n       = grads.nloc(); // nb columns to compute on this PE
    const bool active = a_->active();

    char side;
    char uplo;
    char trans;
    char diag    = 'n';
    double delta = 2.0 * tol;

    if (active)
    {

        std::vector<std::vector<double>> ara(m);
        std::vector<std::vector<double>> aba(m);
        for (int k = 0; k < m; k++)
        {
            ara[k].resize(n);
            aba[k].resize(n);
            for (int i = 0; i < n; i++)
            {
                const int ii = i * nst_ + i + offset_;
                aba[k][i]    = b_[k]->val(ii);
                ara[k][i]    = r_[k]->val(ii);
            }
        }

        std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>*> ra(m);
        std::vector<std::string> names(2 * NDIM);
        names[0] = "matrixRA_0";
        names[1] = "matrixRA_1";
        names[2] = "matrixRA_2";
        names[3] = "matrixRA_3";
#if NDIM > 2
        names[4] = "matrixRA_4";
        names[5] = "matrixRA_5";
#endif
        for (int k = 0; k < m; k++)
        {
            ra[k] = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                names[k], *bcr_, nst_, nst_, nst_, bsize_);
        }

        // build b_ <- a0**-T * b_ * a0**-1
        // a_ upper triangular
        side  = 'r';
        uplo  = 'u';
        trans = 'n';
        for (int k = 0; k < m; k++)
        {
            // b <- b * a0**-1
            b_[k]->trsm(side, uplo, trans, diag, 1., *a_);
        }
        side  = 'l';
        trans = 't';
        for (int k = 0; k < m; k++)
        {
            // b <- a0**-T * b
            b_[k]->trsm(side, uplo, trans, diag, 1., *a_);
        }

        // build r_ <- a0**-T * r_ * a0**-1
        side  = 'r';
        trans = 'n';
        for (int k = 0; k < m; k++)
        {
            // r <- r * a0**-1
            r_[k]->trsm(side, uplo, trans, diag, 1., *a_);
        }
        side  = 'l';
        trans = 't';
        for (int k = 0; k < m; k++)
        {
            // r_ <- a0**-T * r
            r_[k]->trsm(side, uplo, trans, diag, 1., *a_);
        }

        std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>*> ba(m);
        names[0] = "matrixBA_0";
        names[1] = "matrixBA_1";
        names[2] = "matrixBA_2";
        names[3] = "matrixBA_3";
#if NDIM > 2
        names[4] = "matrixBA_4";
        names[5] = "matrixBA_5";
#endif
        for (int k = 0; k < m; k++)
        {
            ba[k] = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                names[k], *bcr_, nst_, nst_, nst_, bsize_);
        }

#ifdef DEBUG
        for (int k = 0; k < m; k++)
        {
            (*MPIdata::sout) << "r_[" << k << "]:" << std::endl;
            r_[k]->print((*MPIdata::sout));
        }
        for (int k = 0; k < m; k++)
        {
            (*MPIdata::sout) << "b_[" << k << "]:" << std::endl;
            b_[k]->print((*MPIdata::sout));
        }
#endif

        assert(tol > 0.0);
        int isweep     = 0;
        double tspread = 0.;

        double oldtspread;

        double trial_step = getDt();

#ifdef DEBUG
        (*MPIdata::sout) << "a initial:" << std::endl;
        a_->print((*MPIdata::sout));
#endif

        getNorms(alphan, inv_alphan);

        for (int k = 0; k < m; k++)
        {
            ra[k]->gemm('n', 'n', 1., *r_[k], *a_, 0.);
            ba[k]->gemm('n', 'n', 1., *b_[k], *a_, 0.);
        }

        // evaluate functionals
        tspread = getFunctionalValue(ara, aba, inv_alphan);

        assert(tspread >= 0.);
        (*MPIdata::sout) << " initial functional value: " << tspread
                         << std::endl;

        // iterations
        while (isweep < maxsweep && delta > tol)
        {

            // compute gradient of spread functional
            grads.clear();
            for (int k = 0; k < m; k++)
            {

                // loop over columns
                for (int i = 0; i < n; i++)
                {
                    grads.axpyColumn(i, alphan[i], *ba[k]);

                    double alpha = ara[k][i] * inv_alphan[i];
                    double fac   = -2. * alpha * alphan[i];
                    grads.axpyColumn(i, fac, *ra[k]);

                    fac = (-1. * aba[k][i] + 2. * alpha * alpha * alphan[i]);
                    grads.axpyColumn(i, fac, *a_);
                }
            }

            for (int i = 0; i < n; i++)
            {
                double fac = 2. * inv_alphan[i];
                grads.scalColumn(i, fac);
            }

            // update the nonorthogonal transformation a_
            const double fac = -1. * trial_step;
            a_->axpy(fac, grads);

            double deriv = 0.;
            for (int i = 0; i < n; i++)
                deriv += grads.dotColumns(i, i) * inv_alphan[i];
            if (n < nst_)
            {
                double recvbf;
                MPI_Allreduce(&deriv, &recvbf, 1, MPI_DOUBLE, MPI_SUM, comm_);
                deriv = recvbf;
            }
            deriv *= fac;

            // update norms
            getNorms(alphan, inv_alphan);

            // update ara and aba
            for (int k = 0; k < m; k++)
            {
                ra[k]->gemm('n', 'n', 1., *r_[k], *a_, 0.);
                for (int i = 0; i < n; i++)
                    ara[k][i] = a_->dotColumns(i, *ra[k], i);

                ba[k]->gemm('n', 'n', 1., *b_[k], *a_, 0.);
                for (int i = 0; i < n; i++)
                    aba[k][i] = a_->dotColumns(i, *ba[k], i);
            }

            // evaluate functionals
            oldtspread = tspread;
            tspread    = getFunctionalValue(ara, aba, inv_alphan);
            //(*MPIdata::sout) << " average spread (trial step):   " <<
            // sqrt(tspread/nst_) << std::endl;

            delta = tspread - oldtspread;

            const double c2 = (delta - deriv);
            const double t  = -deriv / (2. * c2);
            //(*MPIdata::sout) << " new trial t:    " << -1.*t*fac << std::endl;

            double beta = (t - 1.) * fac;
            a_->axpy(beta, grads);

            // update norms
            getNorms(alphan, inv_alphan);

            // update ara and aba
            for (int k = 0; k < m; k++)
            {
                ra[k]->gemm('n', 'n', 1., *r_[k], *a_, 0.);
                for (int i = 0; i < n; i++)
                    ara[k][i] = a_->dotColumns(i, *ra[k], i);

                ba[k]->gemm('n', 'n', 1., *b_[k], *a_, 0.);
                for (int i = 0; i < n; i++)
                    aba[k][i] = a_->dotColumns(i, *ba[k], i);
            }

            // evaluate functional again
            tspread = getFunctionalValue(ara, aba, inv_alphan);

            delta = tspread - oldtspread;

            delta = fabs(delta);
            assert(tspread >= 0.);
            assert(nst_ > 0);
            (*MPIdata::sout) << " functional value: " << tspread << std::endl;

            isweep++;
        } // isweep

        (*MPIdata::sout) << std::setprecision(12);
        (*MPIdata::sout) << isweep << " sweeps in nolmo" << std::endl;
        (*MPIdata::sout) << " tspread:   " << tspread << std::endl;
        (*MPIdata::sout) << " delta: " << delta << std::endl;

        // update r_ and b_ as ara and aba
        for (int k = 0; k < m; k++)
        {
            ra[k]->gemm('n', 'n', 1., *r_[k], *a_, 0.);
            r_[k]->gemm('t', 'n', 1., *a_, *ra[k], 0.);

            ba[k]->gemm('n', 'n', 1., *b_[k], *a_, 0.);
            b_[k]->gemm('t', 'n', 1., *a_, *ba[k], 0.);
        }

#ifdef DEBUG
        (*MPIdata::sout) << "a_ final:" << std::endl;
        a_->print((*MPIdata::sout));
#endif

        for (int k = 0; k < m; k++)
        {
            delete ra[k];
            delete ba[k];
        }
    }

    // get transformation matrix to apply to orbitals
    // to get nonorthogonal localized orbitals:
    // mat_ = a0**-1 * a_ with a0 = ls_**T
    side  = 'l';
    trans = 'n';
    uplo  = 'u';

    gatherTransformMat();

    Ttrsm(side, uplo, trans, diag, nst_, nst_, one, &lst_[0], nst_, &mat_[0],
        nst_);

    return delta;
}

////////////////////////////////////////////////////////////////////////////////

double NOLMOTransform::nolmo_fixedCenters(const int maxsweep, const double tol)
{
    // The input matrices are given in r_[k] and b_[k], k = 0,..,m-1
    // the initial non-orthogonal transformation is stored in a_
    // the non-orthogonal transformation is stored in mat_[i+nst_*j]
    const int m = 2 * NDIM;
    double one  = 1.;

    std::vector<double> alphan(lnst_);
    std::vector<double> inv_alphan(lnst_);

    dist_matrix::DistMatrix<DISTMATDTYPE> grads(
        "grad", *bcr_, nst_, nst_, nst_, bsize_);
    const int n       = grads.nloc(); // nb columns to compute on this PE
    const bool active = a_->active();

    char side;
    char uplo;
    char trans;
    char diag    = 'n';
    double delta = 2.0 * tol;

    if (active)
    {
        std::vector<std::vector<double>> ara(m);
        std::vector<std::vector<double>> aba(m);
        for (int k = 0; k < m; k++)
        {
            ara[k].resize(n);
            aba[k].resize(n);
            for (int i = 0; i < n; i++)
            {
                const int ii = i * nst_ + i + offset_;
                aba[k][i]    = b_[k]->val(ii);
                ara[k][i]    = r_[k]->val(ii);
            }
        }

        std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>*> ra(m);
        std::vector<std::string> names(2 * NDIM);
        names[0] = "matrixRA_0";
        names[1] = "matrixRA_1";
        names[2] = "matrixRA_2";
        names[3] = "matrixRA_3";
#if NDIM > 2
        names[4] = "matrixRA_4";
        names[5] = "matrixRA_5";
#endif
        for (int k = 0; k < m; k++)
        {
            ra[k] = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                names[k], *bcr_, nst_, nst_, nst_, bsize_);
        }

        // build b_ <- a0**-T * b_ * a0**-1
        // a_ upper triangular
        side  = 'r';
        uplo  = 'u';
        trans = 'n';
        for (int k = 0; k < m; k++)
        {
            // b <- b * a0**-1
            b_[k]->trsm(side, uplo, trans, diag, 1., *a_);
        }
        side  = 'l';
        trans = 't';
        for (int k = 0; k < m; k++)
        {
            // b <- a0**-T * b
            b_[k]->trsm(side, uplo, trans, diag, 1., *a_);
        }

        // build r_ <- a0**-T * r_ * a0**-1
        side  = 'r';
        trans = 'n';
        for (int k = 0; k < m; k++)
        {
            // r <- r * a0**-1
            r_[k]->trsm(side, uplo, trans, diag, 1., *a_);
        }
        side  = 'l';
        trans = 't';
        for (int k = 0; k < m; k++)
        {
            // r_ <- a0**-T * r
            r_[k]->trsm(side, uplo, trans, diag, 1., *a_);
        }

        std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>*> ba(m);
        names[0] = "matrixBA_0";
        names[1] = "matrixBA_1";
        names[2] = "matrixBA_2";
        names[3] = "matrixBA_3";
#if NDIM > 2
        names[4] = "matrixBA_4";
        names[5] = "matrixBA_5";
#endif
        for (int k = 0; k < m; k++)
        {
            ba[k] = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                names[k], *bcr_, nst_, nst_, nst_, bsize_);
        }

#ifdef DEBUG
        (*MPIdata::sout) << "a:" << std::endl;
        a_->print((*MPIdata::sout));
        for (int k = 0; k < m; k++)
        {
            (*MPIdata::sout) << "r_[" << k << "]:" << std::endl;
            r_[k]->print((*MPIdata::sout));
        }
        for (int k = 0; k < m; k++)
        {
            (*MPIdata::sout) << "b_[" << k << "]:" << std::endl;
            b_[k]->print((*MPIdata::sout));
        }
#endif

        assert(tol > 0.0);

        double trial_step = getDt();

#ifdef DEBUG
        (*MPIdata::sout) << "a initial:" << std::endl;
        a_->print((*MPIdata::sout));
#endif

        getNorms(alphan, inv_alphan);

        for (int k = 0; k < m; k++)
        {
            ra[k]->gemm('n', 'n', 1., *r_[k], *a_, 0.);
            ba[k]->gemm('n', 'n', 1., *b_[k], *a_, 0.);
        }

        // evaluate functionals
        double tspread = getFunctionalValueFixedCenters(ara, aba, inv_alphan);
        double oldtspread = 0.;

        assert(tspread >= 0.);
        if (bcr_->onpe0())
            (*MPIdata::sout)
                << " initial functional value with fixed centers: " << tspread
                << std::endl;

        // iterations
        int isweep = 0;
        while (isweep < maxsweep && delta > tol)
        {

            // compute gradient of spread functional
            grads.clear();
            for (int k = 0; k < m; k++)
            {
                // loop over columns
                for (int i = 0; i < n; i++)
                {
                    grads.axpyColumn(i, alphan[i], *ba[k]);

                    double alpha = ara0_[k][i];
                    double fac   = -2. * alpha * alphan[i];
                    grads.axpyColumn(i, fac, *ra[k]);

                    fac = (-1. * aba[k][i] + 2. * ara[k][i] * alpha);
                    grads.axpyColumn(i, fac, *a_);
                }
            }

            for (int i = 0; i < n; i++)
            {
                double fac = 2. * inv_alphan[i];
                grads.scalColumn(i, fac);
            }

#ifdef DEBUG
            ra[0]->gemm('t', 'n', 1., grads, *a_, 0.);
            (*MPIdata::sout) << "Test orthogonality of gradient:" << std::endl;
            ra[0]->print((*MPIdata::sout));
#endif

            // update the nonorthogonal transformation a_
            const double fac = -1. * trial_step;
            a_->axpy(fac, grads);

            double deriv = 0.;
            for (int i = 0; i < n; i++)
                deriv += grads.dotColumns(i, i) * inv_alphan[i];
            if (n < nst_)
            {
                double recvbf;
                MPI_Allreduce(&deriv, &recvbf, 1, MPI_DOUBLE, MPI_SUM, comm_);
                deriv = recvbf;
            }
            deriv *= fac;

            // update norms
            getNorms(alphan, inv_alphan);

            // update ara and aba
            for (int k = 0; k < m; k++)
            {
                ra[k]->gemm('n', 'n', 1., *r_[k], *a_, 0.);
                for (int i = 0; i < n; i++)
                    ara[k][i] = a_->dotColumns(i, *ra[k], i);

                ba[k]->gemm('n', 'n', 1., *b_[k], *a_, 0.);
                for (int i = 0; i < n; i++)
                    aba[k][i] = a_->dotColumns(i, *ba[k], i);
            }

            // evaluate functionals
            oldtspread = tspread;
            tspread    = getFunctionalValueFixedCenters(ara, aba, inv_alphan);
            //(*MPIdata::sout) << " average spread (trial step):   " <<
            // sqrt(tspread/nst_) << std::endl;

            delta = tspread - oldtspread;

            const double c2 = (delta - deriv);
            const double t  = -deriv / (2. * c2);
            //(*MPIdata::sout) << " new trial t: " << -1.*t*fac << std::endl;

            double beta = (t - 1.) * fac;
            a_->axpy(beta, grads);

            // update norms
            getNorms(alphan, inv_alphan);

            // update ara and aba
            for (int k = 0; k < m; k++)
            {
                ra[k]->gemm('n', 'n', 1., *r_[k], *a_, 0.);
                for (int i = 0; i < n; i++)
                    ara[k][i] = a_->dotColumns(i, *ra[k], i);

                ba[k]->gemm('n', 'n', 1., *b_[k], *a_, 0.);
                for (int i = 0; i < n; i++)
                    aba[k][i] = a_->dotColumns(i, *ba[k], i);
            }

            // evaluate functional again
            tspread = getFunctionalValueFixedCenters(ara, aba, inv_alphan);

            delta = tspread - oldtspread;

            delta = fabs(delta);
            assert(tspread >= 0.);
            assert(nst_ > 0);
            if (bcr_->onpe0())
                (*MPIdata::sout)
                    << " functional value with fixed centers: " << tspread
                    << std::endl;

            isweep++;
        } // isweep

        if (bcr_->onpe0())
        {
            (*MPIdata::sout) << std::setprecision(12);
            (*MPIdata::sout)
                << isweep << " sweeps in nolmo_fixedCenters" << std::endl;
            (*MPIdata::sout) << " tspread:   " << tspread << std::endl;
            (*MPIdata::sout) << " delta: " << delta << std::endl;
        }

        // update r_ and b_ as ara and aba
        for (int k = 0; k < m; k++)
        {
            ra[k]->gemm('n', 'n', 1., *r_[k], *a_, 0.);
            r_[k]->gemm('t', 'n', 1., *a_, *ra[k], 0.);

            ba[k]->gemm('n', 'n', 1., *b_[k], *a_, 0.);
            b_[k]->gemm('t', 'n', 1., *a_, *ba[k], 0.);
        }

#ifdef DEBUG
        (*MPIdata::sout) << "a_ final:" << std::endl;
        a_->print((*MPIdata::sout));
#endif

        for (int k = 0; k < m; k++)
        {
            delete ra[k];
            delete ba[k];
        }

    } // active

    // get transformation matrix to apply to orbitals
    // to get nonorthogonal localized orbitals:
    // mat_ = a0**-1 * a_ with a0 = ls_**T
    side  = 'l';
    trans = 'n';
    uplo  = 'u';
    // a_->matgather(&mat_[0],nst_);

    gatherTransformMat();

#if 0
  if( nst_<20 ){
    (*MPIdata::sout)<<"NOLMO Transformation matrix (before):"<<std::endl;
    (*MPIdata::sout)<<std::setprecision(2)<<scientific;
    for ( int i = 0; i < nst_; i++ ){
      for ( int j = 0; j < nst_; j++ ){
        (*MPIdata::sout)<<mat_[i+nst_*j]<<"\t";
      }
      (*MPIdata::sout)<<std::endl;
    }
  }
#endif
    Ttrsm(side, uplo, trans, diag, nst_, nst_, one, &lst_[0], nst_, &mat_[0],
        nst_);

#if 0
  if( nst_<20 ){
    (*MPIdata::sout)<<"NOLMO Transformation matrix:"<<std::endl;
    (*MPIdata::sout)<<std::setprecision(2)<<scientific;
    for ( int i = 0; i < nst_; i++ ){
      for ( int j = 0; j < nst_; j++ ){
        (*MPIdata::sout)<<mat_[i+nst_*j]<<"\t";
      }
      (*MPIdata::sout)<<std::endl;
    }
  }
#endif

    return delta;
}

void NOLMOTransform::gatherTransformMat()
{
    const int nloc   = (a_->active()) ? a_->nloc() : 0;
    const int nbccol = a_->npcol();
    const int npes   = bcr_->nprocs();

    std::vector<int> recvcounts(npes, 0);
    for (int i = 0; i < nbccol; i++)
    {
        if ((i + 1) * bsize_ < nst_)
            recvcounts[i] = bsize_ * nst_;
        else
            recvcounts[i] = std::max(0, nst_ * (nst_ - i * bsize_));
    }
    for (unsigned int i = 0; i < recvcounts.size(); i++)
    {
        assert(recvcounts[i] <= bsize_ * nst_);
        assert(recvcounts[i] >= 0);
    }
    std::vector<int> displs(npes, 0);
    for (int i = 0; i < nbccol; i++)
        displs[i] = i * bsize_ * nst_;
    // The block of data sent from the jth process is received by every
    // process and placed in the jth block of the buffer recvbuf.
    std::vector<DISTMATDTYPE> sendbuf(nst_ * nloc);
    if (a_->active()) a_->copyDataToVector(sendbuf);
    DISTMATDTYPE* recvbuf = &mat_[0];
    assert(recvbuf != nullptr);
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();
    MPI_Allgatherv(&sendbuf[0], nst_ * nloc, MPI_DOUBLE, recvbuf,
        &recvcounts[0], &displs[0], MPI_DOUBLE, comm);
}

double NOLMOTransform::getFunctionalValueFixedCenters(
    const std::vector<std::vector<double>>& ara,
    const std::vector<std::vector<double>>& aba,
    const std::vector<double>& inv_alphan)
{
    assert(ara0_[0].size() == ara[0].size());
    assert(ara0_[0].size() == aba[0].size());
    assert(ara0_[0].size() == inv_alphan.size());

    double tspread = 0.0;
    int n          = (int)ara0_[0].size();
    for (int i = 0; i < n; i++)
    {
        double spread2 = 0.0;
        for (int k = 0; k < 2 * NDIM; k++)
        {
            const double vtmp  = ara[k][i] * inv_alphan[i];
            const double vtmp0 = ara0_[k][i];
            spread2 += ((aba[k][i] * inv_alphan[i]) + vtmp0 * vtmp0
                        - 2. * vtmp0 * vtmp);
        }
        //(*MPIdata::sout) << " spread2["<<i<<"]:   " << spread2 << std::endl;
        tspread += spread2;
    }

    if (bcr_->npcol() > 1)
    {
        //(*MPIdata::sout) << "tspread="<<tspread<<std::endl;
        double recvbf;
        MPI_Allreduce(&tspread, &recvbf, 1, MPI_DOUBLE, MPI_SUM, comm_);
        tspread = recvbf;
    }
    return tspread;
}

double NOLMOTransform::getFunctionalValue(
    const std::vector<std::vector<double>>& ara,
    const std::vector<std::vector<double>>& aba,
    const std::vector<double>& inv_alphan)
{
    assert(ara0_[0].size() == ara[0].size());
    assert(ara0_[0].size() == aba[0].size());
    assert(ara0_[0].size() == inv_alphan.size());

    double tspread = 0.0;
    int n          = (int)ara0_[0].size();
    for (int i = 0; i < n; i++)
    {
        double spread2 = 0.0;
        for (int k = 0; k < 2 * NDIM; k++)
        {
            const double vtmp = ara[k][i] * inv_alphan[i];
            spread2 += ((aba[k][i] * inv_alphan[i]) - vtmp * vtmp);
        }
        //(*MPIdata::sout) << " spread2["<<i<<"]:   " << spread2 << std::endl;
        tspread += spread2;
    }

    if (bcr_->npcol() > 1)
    {
        double recvbf;
        MPI_Allreduce(&tspread, &recvbf, 1, MPI_DOUBLE, MPI_SUM, comm_);
        tspread = recvbf;
    }
    return tspread;
}

void NOLMOTransform::getNorms(
    std::vector<double>& alphan, std::vector<double>& inv_alphan)
{
    int n = (int)inv_alphan.size();

    // sweep: loop over columns to get norms of vectors
    for (int p = 0; p < n; p++)
    {
        alphan[p] = a_->dotColumns(p, p);
        assert(alphan[p] > 0.);
        inv_alphan[p] = 1. / alphan[p];
    }
}

double NOLMOTransform::getDt()
{
    double max_cell = cell_[0];
    for (int i = 1; i < NDIM; i++)
        if (cell_[i] > max_cell) max_cell = cell_[i];
    double trial_step = 2. * M_PI / max_cell;
    trial_step *= trial_step;
    trial_step /= 12.;

    //(*MPIdata::sout)<<"trial_step="<<trial_step<<std::endl;

    return trial_step;
}
