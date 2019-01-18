// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
//#define DEBUG 1
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Control.h"
#include "MGmol_MPI.h"
#include "MLWFTransform.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef ADD_
#define drot drot_
#define srot srot_
#endif

extern "C"
{
    void drot(int*, double*, int*, double*, int*, double*, double*);
    void srot(int*, float*, int*, float*, int*, float*, float*);
}

using namespace std;

//#if 0
////////////////////////////////////////////////////////////////////////////////
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

#ifdef USE_MPI
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
#endif

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

    double delta = jade(maxsweep, tol);

    if (bcr_->onpe0())
        if (delta > tol)
            (*MPIdata::sout)
                << " MLWFTransform: decrease was " << setprecision(12) << delta
                << " after " << maxsweep << " iterations" << endl;
}

////////////////////////////////////////////////////////////////////////////////

double MLWFTransform::jade(int maxsweep, double tol)
{
    assert(tol > 0.0);

    const int m = 2 * NDIM;

    const int npecol = a_->npcol();
#ifdef SCALAPACK
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int mype        = mmpi.mypeSpin();

    MPI_Status status[8];
    MPI_Request req[8];

    assert(a_->nprow() == 1);
    assert(a_->npcol() <= mmpi.size());
#else
    const int mype = 0;
#endif

    // joint approximate diagonalization of a set
    // of m real symmetric nst_*n matrices
    // The input matrices are given in r_[k], k = 0,..,m-1
    // the element a(i,j) of matrix a_k is in r_[k][i+nst_*j];
    // the orthogonal transformation is returned in mat_[i+nst_*j]

    double sum    = 0.0; // sum of squares of diagonal elements
    double oldsum = 0.0;
    double delta  = 2.0 * tol;

    int isweep = 0;
    int n_loc  = bsize_;

    const int n_loc_2 = bsize_ / 2;
    const int nn      = npecol * bsize_;
    //(*MPIdata::sout)<<"n_loc="<<n_loc<<endl;
    //(*MPIdata::sout)<<"n_loc_2="<<n_loc_2<<endl;
    //(*MPIdata::sout)<<"mype="<<mype<<", offset_="<<offset_<<endl;

    // local storage for matrix u
    vector<double> u_loc(nst_ * bsize_, 0.);

    // pointers to the actual index of the columns
    // actual[j]>=nst_ correspond to dummy columns
    vector<int> actual(bsize_, nst_);
    if (mype < npecol)
    {

        for (int j = 0; j < bsize_; j++)
            actual[j] = j + offset_;

        vector<double> sbuf1((m + 1) * nst_);
        vector<double> sbuf2((m + 1) * nst_);
        vector<double> rbuf1((m + 1) * nst_);
        vector<double> rbuf2((m + 1) * nst_);

        const double eps = 1.e-14;

        vector<double> u_small(bsize_);

        // local storage for matrices r_[k]
        vector<vector<DISTMATDTYPE>> r_loc;
        r_loc.resize(m);
        for (int k = 0; k < m; k++)
            r_loc[k].resize(nst_ * bsize_);

        // initialize the orthogonal transformation with the identity operator
        for (int i = 0; i < bsize_; i++)
            if (actual[i] < nst_) u_loc[actual[i] + nst_ * i] = 1.0;

        // initialize r_loc:
        for (int k = 0; k < m; k++)
            memset(&r_loc[k][0], 0, nst_ * bsize_ * sizeof(DISTMATDTYPE));
        for (int k = 0; k < m; k++)
            for (int j = 0; j < bsize_; j++)
            {
                const int nj = nst_ * j;
                if (actual[j] < nst_)
                    for (int i = 0; i < nst_; i++)
                    {
                        r_loc[k][nj + i] = r_[k]->val(nj + i);
                    }
            }

        // indexing current local columns
        vector<int> column(bsize_);
        for (int i = 0; i < bsize_; i++)
            column[i] = i;

#if DEBUG
        if (bcr_->onpe0()) (*MPIdata::sout) << "r_loc" << endl;
        if (bsize_ < 8)
            for (int k = 0; k < m; k++)
                for (int j = 0; j < bsize_; j++)
                {
                    for (int i = 0; i < nst_; i++)
                        (*MPIdata::sout) << r_loc[k][i + nst_ * j] << "\t";
                    (*MPIdata::sout) << endl;
                }
#endif

        while (isweep < maxsweep && delta > tol)
        {
            // do nn-1 rounds
            for (int round = 0; round < nn - 1; round++)
            {
#if DEBUG
                // if(bcr_->onpe0())
                (*MPIdata::sout) << "JADE: round " << round << endl;
#endif

                vector<int> actualv(bsize_);

                // loop over local pairs
                for (int ip = 0; ip < n_loc_2; ip++)
                {
                    const int ip2 = 2 * ip;
                    const int iq  = bsize_ - 1 - ip; // iq paired with ip
                    const int cp  = column[ip];
                    const int cq  = column[iq];
                    const int rp  = actual[cp];
                    const int rq  = actual[cq];

                    // compute the rotation minimizing sum_k a_pp^2 + a_qq^2

                    if (rp < nst_ && rq < nst_)
                    // if( iu_[rp+nst_*rq] )
                    {

#if DEBUG
                        (*MPIdata::sout) << "compute actual pair: " << rp
                                         << ", " << rq << endl;
#endif
                        double g11 = 0.0, g12 = 0.0, g22 = 0.0;

                        for (int k = 0; k < m; k++)
                        {
                            double app = 1.;
                            double aqq = 0.;
                            double apq = 0.;
                            if (rp < nst_ && rq < nst_)
                            {
                                apq = r_loc[k][rp + nst_ * cq];
                                app = r_loc[k][rp + nst_ * cp];
                                aqq = r_loc[k][rq + nst_ * cq];
                                if (!(fabs(r_loc[k][rp + nst_ * cq]
                                           - r_loc[k][rq + nst_ * cp])
                                        < 0.00001))
                                {
                                    (*MPIdata::sout) << "mype=" << mype
                                                     << ", k=" << k << endl;
                                    (*MPIdata::sout)
                                        << "mype=" << mype
                                        << ", r_loc: rp=" << rp << ", cq=" << cq
                                        << ", " << r_loc[k][rp + nst_ * cq]
                                        << endl;
                                    (*MPIdata::sout)
                                        << "mype=" << mype
                                        << ", r_loc: rq=" << rq << ", cp=" << cp
                                        << ", " << r_loc[k][rq + nst_ * cp]
                                        << endl;
                                }
                                assert(fabs(r_loc[k][rp + nst_ * cq]
                                            - r_loc[k][rq + nst_ * cp])
                                       < 0.00001);
                            }

                            const double h1 = app - aqq;
                            const double h2 = 2.0 * apq;

                            g11 += h1 * h1;
                            g12 += h1 * h2;
                            g22 += h2 * h2;
                        }

                        // the matrix g is real, symmetric

                        // compute Jacobi rotation diagonalizing g

                        double c = 1.0, s = 0.0, e1 = g11, e2 = g22;
                        if (g12 * g12 > eps * eps * fabs(g11 * g22))
                        {
                            double tau = 0.5 * (g22 - g11) / g12;
                            double t
                                = 1.0 / (fabs(tau) + sqrt(1.0 + tau * tau));
                            if (tau < 0.0) t *= -1.0;
                            c = 1.0 / sqrt(1.0 + t * t);
                            s = t * c;

                            // eigenvalues
                            e1 -= t * g12;
                            e2 += t * g12;
                        }

                        // components of eigenvector associated with the largest
                        // eigenvalue
                        double x, y;
                        if (e1 > e2)
                        {
                            x = c;
                            y = -s;
                        }
                        else
                        {
                            x = s;
                            y = c;

                            // for "inner" rotations
                            if (x < 0.)
                            {
                                x = -x;
                                y = -y;
                            }
                        }

                        // compute Jacobi rotation R(p,q)
                        c = sqrt(0.5 * (x + 1.0));
                        s = y / sqrt(2.0 * (x + 1.0));

                        const int ncp = nst_ * cp;
                        const int ncq = nst_ * cq;

                        // update the orthogonal transformation u with the 2*2
                        // rotation R(p,q)
                        //
                        //           |  c   s |
                        //  R(p,q) = |        |
                        //           | -s   c |

                        // U := U * R(p,q)^T
#if DEBUG
                        (*MPIdata::sout) << "c=" << c << ", s=" << s << endl;
#endif
                        int ione = 1;
                        drot(&nst_, &u_loc[ncp], &ione, &u_loc[ncq], &ione, &c,
                            &s);

                        // update the matrices r_loc[k]
                        for (int k = 0; k < m; k++)
                        {
                            // A := A * R(p,q)^T
                            int ione = 1;
                            drot(&nst_, &r_loc[k][ncp], &ione, &r_loc[k][ncq],
                                &ione, &c, &s);
                        }

                        u_small[ip2]     = c;
                        u_small[ip2 + 1] = s;
                    }
                    actualv[ip2]     = rp;
                    actualv[ip2 + 1] = rq;
                } // loop over local pairs

                vector<double> v_small(n_loc);
                vector<int> remote_actualv(n_loc);
                for (int pe = 0; pe < npecol; pe++)
                {

#ifdef SCALAPACK
                    if (mype == pe)
                    {
                        v_small        = u_small;
                        remote_actualv = actualv;
                    }

                    MPI_Bcast(&v_small[0], n_loc, MPI_DOUBLE, pe, comm_);
                    MPI_Bcast(&remote_actualv[0], n_loc, MPI_INT, pe, comm_);
#else
                    v_small        = u_small;
                    remote_actualv = actualv;
#endif
                    // Rotate rows rp and rq
                    for (int ip = 0; ip < n_loc_2; ip++)
                    {
                        const int ip2 = 2 * ip;

                        // rows to rotate
                        const int rowp = remote_actualv[ip2];
                        const int rowq = remote_actualv[ip2 + 1];

                        // rotation coefficients
                        double c = v_small[ip2];
                        double s = v_small[ip2 + 1];

                        // update the rows of the matrices r_loc[k]
                        if (rowp < nst_ && rowq < nst_)
                        // if( iu_[rowp+nst_*rowq] )
                        {
                            for (int k = 0; k < m; k++)
                            {
                                // A := R(p,q) * A
                                drot(&n_loc, &r_loc[k][rowp], &nst_,
                                    &r_loc[k][rowq], &nst_, &c, &s);
                            }
                        }
                    } // ip

                } // pe

#ifdef SCALAPACK
                for (int i = 0; i < 8; i++)
                    req[i] = MPI_REQUEST_NULL;

                int destright = (mype + 1) % npecol;
                int destleft  = (mype - 1) % npecol;

                // save data to exchange in buffers
                assert(column[(n_loc_2 - 1)] < n_loc);
                assert(column[(n_loc - 1)] < n_loc);

                // copy column n_loc_2-1 into buffer to transfer right
                if (mype != npecol - 1)
                {
                    for (int k = 0; k < m; k++)
                        memcpy(&sbuf1[k * nst_],
                            &r_loc[k][column[(n_loc_2 - 1)] * nst_],
                            nst_ * sizeof(DISTMATDTYPE));
                    memcpy(&sbuf1[m * nst_],
                        &u_loc[column[(n_loc_2 - 1)] * nst_],
                        nst_ * sizeof(double));
                }
                // copy column n_loc-1 into buffer to transfer left
                if (mype != 0)
                {
                    for (int k = 0; k < m; k++)
                        memcpy(&sbuf2[k * nst_],
                            &r_loc[k][column[(n_loc - 1)] * nst_],
                            nst_ * sizeof(double));
                    memcpy(&sbuf2[m * nst_], &u_loc[column[(n_loc - 1)] * nst_],
                        nst_ * sizeof(double));
                }
                int scolr = actual[column[(n_loc_2 - 1)]];
                int scoll = actual[column[(n_loc - 1)]];

                // Exchange columns of matrices
                if (mype != npecol - 1)
                {
                    MPI_Irecv(&rbuf2[0], (m + 1) * nst_, MPI_DOUBLE, destright,
                        destright, comm_, req);
                    MPI_Irecv(&actual[column[(n_loc_2 - 1)]], 1, MPI_INT,
                        destright, destright + npecol, comm_, req + 4);
                }

                if (mype != 0)
                {
                    MPI_Irecv(&rbuf1[0], (m + 1) * nst_, MPI_DOUBLE, destleft,
                        destleft, comm_, req + 1);
                    MPI_Irecv(&actual[column[(n_loc - 1)]], 1, MPI_INT,
                        destleft, destleft + npecol, comm_, req + 5);
                }

                if (mype != 0)
                {
                    MPI_Isend(&sbuf2[0], (m + 1) * nst_, MPI_DOUBLE, destleft,
                        mype, comm_, req + 2);
                    MPI_Isend(&scoll, 1, MPI_INT, destleft, mype + npecol,
                        comm_, req + 6);
                    //(*MPIdata::sout)<<"mype="<<mype<<", send col left
                    //"<<scoll<<endl;
                }

                if (mype != npecol - 1)
                {
                    MPI_Isend(&sbuf1[0], (m + 1) * nst_, MPI_DOUBLE, destright,
                        mype, comm_, req + 3);
                    MPI_Isend(&scolr, 1, MPI_INT, destright, mype + npecol,
                        comm_, req + 7);
                    //(*MPIdata::sout)<<"mype="<<mype<<", send col right
                    //"<<scolr<<endl;
                }

                MPI_Wait(req, status);
                MPI_Wait(req + 1, status + 1);
                MPI_Wait(req + 2, status + 2);
                MPI_Wait(req + 3, status + 3);
                MPI_Wait(req + 4, status + 4);
                MPI_Wait(req + 5, status + 5);
                MPI_Wait(req + 6, status + 6);
                MPI_Wait(req + 7, status + 7);

                for (int k = 0; k < m; k++)
                {
                    if (mype != npecol - 1)
                        memcpy(&r_loc[k][column[(n_loc_2 - 1)] * nst_],
                            &rbuf2[k * nst_], nst_ * sizeof(double));
                    if (mype != 0)
                        memcpy(&r_loc[k][column[(n_loc - 1)] * nst_],
                            &rbuf1[k * nst_], nst_ * sizeof(double));
                }
                if (mype != npecol - 1)
                    memcpy(&u_loc[column[(n_loc_2 - 1)] * nst_],
                        &rbuf2[m * nst_], nst_ * sizeof(double));
                if (mype != 0)
                    memcpy(&u_loc[column[(n_loc - 1)] * nst_], &rbuf1[m * nst_],
                        nst_ * sizeof(double));
#endif

                // rotate pointers for local columns
                vector<int>::iterator bp = column.begin();
                if (mype == 0) bp++;
                vector<int>::iterator ep = column.end();
                ep--;
                const int elm = (*ep); // value of last element
                assert(elm < n_loc);
                column.insert(bp, elm);
                column.pop_back();

            } // end round

            oldsum = sum;
            sum    = 0.0; // sum of squares of diagonal elements

            double sum_off = 0.;
            for (int k = 0; k < m; k++)
            {
                // loop over local columns
                for (int ip = 0; ip < n_loc; ip++)
                {
                    int rp = actual[ip];
                    if (rp < nst_)
                    {
                        sum += r_loc[k][rp + nst_ * ip]
                               * r_loc[k][rp + nst_ * ip];
                        for (int i = 0; i < nst_; i++)
                            sum_off += r_loc[k][i + nst_ * ip]
                                       * r_loc[k][i + nst_ * ip];
                        sum_off -= r_loc[k][rp + nst_ * ip]
                                   * r_loc[k][rp + nst_ * ip];
                    }
                }
            }

#ifdef SCALAPACK
            double in;
            int mpirc;

            in    = sum;
            mpirc = MPI_Allreduce(&in, &sum, 1, MPI_DOUBLE, MPI_SUM, comm_);
            in    = sum_off;
            mpirc = MPI_Allreduce(&in, &sum_off, 1, MPI_DOUBLE, MPI_SUM, comm_);
#endif
#if DEBUG
            (*MPIdata::sout) << " Sum_off(A(k)):  " << sum_off << endl;
            (*MPIdata::sout) << " Sum_diag(A(k)): " << sum << endl;
            (*MPIdata::sout) << " Frobenius:        " << sum_off + sum << endl;
            (*MPIdata::sout) << endl;
#endif

            isweep++;

            // delta = sum - oldsum;
            delta = fabs(sum - oldsum);
#if DEBUG
            if (bcr_->onpe0()) (*MPIdata::sout) << " sum:   " << sum << endl;
            if (bcr_->onpe0()) (*MPIdata::sout) << " delta: " << delta << endl;
#endif

        } // isweep

        // loop over local columns
        for (int ip = 0; ip < n_loc; ip++)
        {
#if DEBUG
            // check columns back onto original PE
            if (!(actual[ip] < (mype + 1) * n_loc || actual[ip] == -1))
            {
                (*MPIdata::sout) << "mype=" << mype << ", ip=" << ip
                                 << ", actual[ip]=" << actual[ip] << endl;
            }
            if (!(actual[ip] >= mype * n_loc || actual[ip] == -1))
            {
                (*MPIdata::sout) << "mype=" << mype << ", ip=" << ip
                                 << ", actual[ip]=" << actual[ip] << endl;
            }
#endif
            assert(actual[ip] < (mype + 1) * n_loc || actual[ip] == -1);
            assert(actual[ip] >= mype * n_loc || actual[ip] == -1);

            if (actual[ip] < nst_)
            {
                // update columns ip of r_
                for (int k = 0; k < m; k++)
                {
                    assert(ip < r_[k]->nb());
                    r_[k]->assignColumn(
                        &r_loc[k][nst_ * ip], actual[ip] - offset_);
                }
            }
        }

    } // if mype<npecol

#if DEBUG
    for (int k = 0; k < m; k++)
    {
        if (bcr_->onpe0()) (*MPIdata::sout) << "r_ for k=" << k << endl;
        r_[k]->print((*MPIdata::sout), 0, 0, 5, 5);
    }
#endif

    // bcast u_loc into mat_
    for (int pe = 0; pe < npecol; pe++)
    {
        // bcast actual
        int root = pe;
        // loop over local columns
        for (int ip = 0; ip < n_loc; ip++)
        {
            int lip = actual[ip];
#ifdef SCALAPACK
            mmpi.bcast(&lip, 1, root);
#endif
            if (lip < nst_)
            {
                if (mype == root)
                {
                    memcpy(&mat_[nst_ * lip], &u_loc[nst_ * ip],
                        nst_ * sizeof(double));
                }
#ifdef SCALAPACK
                mmpi.bcast(&mat_[nst_ * lip], nst_, root);
#endif
            }
        }
    }

    if (bcr_->onpe0())
    {
        (*MPIdata::sout) << setprecision(12);
        (*MPIdata::sout) << isweep << " sweeps in jade parallel" << endl;
        (*MPIdata::sout) << " sum:   " << sum << endl;
        (*MPIdata::sout) << " delta: " << delta << endl;
    }

    return delta;
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
//#endif
