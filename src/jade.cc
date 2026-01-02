// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrix.h"
#endif
#include "MGmol_MPI.h"
#include "MGmol_blas1.h"
#include "fc_mangle.h"

#include <cassert>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <mpi.h>

#define NDIM 3

// joint approximate diagonalization of a set
// of m real symmetric m_loc*n matrices
// The input matrices are given in rmat[k], k = 0,..,m-1
// the element a(i,j) of matrix a_k is in rmat[k][i+m_loc*j];
// the orthogonal transformation is returned in tmat[i+m_loc*j]
#ifdef MGMOL_USE_SCALAPACK
double jade(std::vector<dist_matrix::DistMatrix<double>*>& rmat,
    const int n_loc, const int m_loc, const int offset_, const int maxsweep,
    const double tol, MPI_Comm comm, std::vector<double>& tmat,
    const bool print_flag)
{
    assert(tol > 0.0);

    const int m = 2 * NDIM;

    const int npecol = rmat[0]->npcol();
#ifdef MGMOL_USE_SCALAPACK
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int mype        = mmpi.mypeSpin();

    MPI_Status status[8];
    MPI_Request req[8];

    assert(rmat[0]->nprow() == 1);
    assert(rmat[0]->npcol() <= mmpi.size());
#else
    const int mype = 0;
#endif

    double sum    = 0.0; // sum of squares of diagonal elements
    double oldsum = 0.0;
    double delta  = 2.0 * tol;

    int isweep = 0;
    //    int n_loc  = n_loc;

    const int n_loc_2 = n_loc / 2;
    const int nn      = npecol * n_loc;
    // cout<<"n_loc="<<n_loc<<endl;
    // cout<<"n_loc_2="<<n_loc_2<<endl;
    // cout<<"mype="<<mype<<", offset_="<<offset_<<endl;

    // local storage for matrix u
    std::vector<double> u_loc(m_loc * n_loc, 0.);

    // pointers to the actual index of the columns
    // actual[j]>=m_loc correspond to dummy columns
    std::vector<int> actual(n_loc, m_loc);
    if (mype < npecol)
    {
        for (int j = 0; j < n_loc; j++)
            actual[j] = j + offset_;

        std::vector<double> sbuf1((m + 1) * m_loc);
        std::vector<double> sbuf2((m + 1) * m_loc);
        std::vector<double> rbuf1((m + 1) * m_loc);
        std::vector<double> rbuf2((m + 1) * m_loc);

        const double eps = 1.e-14;

        std::vector<double> u_small(n_loc);

        // local storage for matrices rmat[k]
        std::vector<std::vector<double>> r_loc;
        r_loc.resize(m);
        for (int k = 0; k < m; k++)
            r_loc[k].resize(m_loc * n_loc);

        // initialize the orthogonal transformation with the identity operator
        for (int i = 0; i < n_loc; i++)
            if (actual[i] < m_loc) u_loc[actual[i] + m_loc * i] = 1.0;

        // initialize r_loc:
        for (int k = 0; k < m; k++)
            memset(&r_loc[k][0], 0, m_loc * n_loc * sizeof(double));
        for (int k = 0; k < m; k++)
            for (int j = 0; j < n_loc; j++)
            {
                const int nj = m_loc * j;
                if (actual[j] < m_loc)
                    for (int i = 0; i < m_loc; i++)
                    {
                        r_loc[k][nj + i] = rmat[k]->val(nj + i);
                    }
            }

        // indexing current local columns
        std::vector<int> column(n_loc);
        for (int i = 0; i < n_loc; i++)
            column[i] = i;

#if DEBUG
        if (print_flag) std::cout << "r_loc" << std::endl;
        if (n_loc < 8)
            for (int k = 0; k < m; k++)
                for (int j = 0; j < n_loc; j++)
                {
                    for (int i = 0; i < m_loc; i++)
                        std::cout << r_loc[k][i + m_loc * j] << "\t";
                    std::cout << std::endl;
                }
#endif

        int nrot = n_loc;
        int mrot = m_loc;

        while (isweep < maxsweep && delta > tol)
        {
            // do nn-1 rounds
            for (int round = 0; round < nn - 1; round++)
            {
#if DEBUG
                // if(print_flag)
                std::cout << "JADE: round " << round << std::endl;
#endif

                std::vector<int> actualv(n_loc);

                // loop over local pairs
                for (int ip = 0; ip < n_loc_2; ip++)
                {
                    const int ip2 = 2 * ip;
                    const int iq  = n_loc - 1 - ip; // iq paired with ip
                    const int cp  = column[ip];
                    const int cq  = column[iq];
                    const int rp  = actual[cp];
                    const int rq  = actual[cq];

                    // compute the rotation minimizing sum_k a_pp^2 + a_qq^2

                    if (rp < m_loc && rq < m_loc)
                    // if( iu_[rp+m_loc*rq] )
                    {
#if DEBUG
                        std::cout << "compute actual pair: " << rp << ", " << rq
                                  << std::endl;
#endif
                        double g11 = 0.0, g12 = 0.0, g22 = 0.0;

                        for (int k = 0; k < m; k++)
                        {
                            double app = 1.;
                            double aqq = 0.;
                            double apq = 0.;
                            if (rp < m_loc && rq < m_loc)
                            {
                                apq = r_loc[k][rp + m_loc * cq];
                                app = r_loc[k][rp + m_loc * cp];
                                aqq = r_loc[k][rq + m_loc * cq];
                                if (!(std::abs(r_loc[k][rp + m_loc * cq]
                                               - r_loc[k][rq + m_loc * cp])
                                        < 0.00001))
                                {
                                    std::cout << "mype=" << mype << ", k=" << k
                                              << std::endl;
                                    std::cout << "mype=" << mype
                                              << ", r_loc: rp=" << rp
                                              << ", cq=" << cq << ", "
                                              << r_loc[k][rp + m_loc * cq]
                                              << std::endl;
                                    std::cout << "mype=" << mype
                                              << ", r_loc: rq=" << rq
                                              << ", cp=" << cp << ", "
                                              << r_loc[k][rq + m_loc * cp]
                                              << std::endl;
                                }
                                assert(std::abs(r_loc[k][rp + m_loc * cq]
                                                - r_loc[k][rq + m_loc * cp])
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
                        if (g12 * g12 > eps * eps * std::abs(g11 * g22))
                        {
                            double tau = 0.5 * (g22 - g11) / g12;
                            double t
                                = 1.0 / (std::abs(tau) + sqrt(1.0 + tau * tau));
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

                        const int ncp = m_loc * cp;
                        const int ncq = m_loc * cq;

                        // update the orthogonal transformation u with the 2*2
                        // rotation R(p,q)
                        //
                        //           |  c   s |
                        //  R(p,q) = |        |
                        //           | -s   c |

                        // U := U * R(p,q)^T
#if DEBUG
                        std::cout << "c=" << c << ", s=" << s << std::endl;
#endif
                        int ione = 1;
                        DROT(&mrot, &u_loc[ncp], &ione, &u_loc[ncq], &ione, &c,
                            &s);

                        // update the matrices r_loc[k]
                        for (int k = 0; k < m; k++)
                        {
                            // A := A * R(p,q)^T
                            int ione = 1;
                            DROT(&mrot, &r_loc[k][ncp], &ione, &r_loc[k][ncq],
                                &ione, &c, &s);
                        }

                        u_small[ip2]     = c;
                        u_small[ip2 + 1] = s;
                    }
                    actualv[ip2]     = rp;
                    actualv[ip2 + 1] = rq;
                } // loop over local pairs

                std::vector<double> v_small(n_loc);
                std::vector<int> remote_actualv(n_loc);
                for (int pe = 0; pe < npecol; pe++)
                {
#ifdef MGMOL_USE_SCALAPACK
                    if (mype == pe)
                    {
                        v_small        = u_small;
                        remote_actualv = actualv;
                    }

                    MPI_Bcast(&v_small[0], n_loc, MPI_DOUBLE, pe, comm);
                    MPI_Bcast(&remote_actualv[0], n_loc, MPI_INT, pe, comm);
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
                        if (rowp < m_loc && rowq < m_loc)
                        // if( iu_[rowp+m_loc*rowq] )
                        {
                            for (int k = 0; k < m; k++)
                            {
                                // A := R(p,q) * A
                                DROT(&nrot, &r_loc[k][rowp], &mrot,
                                    &r_loc[k][rowq], &mrot, &c, &s);
                            }
                        }
                    } // ip

                } // pe

#ifdef MGMOL_USE_SCALAPACK
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
                        memcpy(&sbuf1[k * m_loc],
                            &r_loc[k][column[(n_loc_2 - 1)] * m_loc],
                            m_loc * sizeof(double));
                    memcpy(&sbuf1[m * m_loc],
                        &u_loc[column[(n_loc_2 - 1)] * m_loc],
                        m_loc * sizeof(double));
                }
                // copy column n_loc-1 into buffer to transfer left
                if (mype != 0)
                {
                    for (int k = 0; k < m; k++)
                        memcpy(&sbuf2[k * m_loc],
                            &r_loc[k][column[(n_loc - 1)] * m_loc],
                            m_loc * sizeof(double));
                    memcpy(&sbuf2[m * m_loc],
                        &u_loc[column[(n_loc - 1)] * m_loc],
                        m_loc * sizeof(double));
                }
                int scolr = actual[column[(n_loc_2 - 1)]];
                int scoll = actual[column[(n_loc - 1)]];

                // Exchange columns of matrices
                if (mype != npecol - 1)
                {
                    MPI_Irecv(&rbuf2[0], (m + 1) * m_loc, MPI_DOUBLE, destright,
                        destright, comm, req);
                    MPI_Irecv(&actual[column[(n_loc_2 - 1)]], 1, MPI_INT,
                        destright, destright + npecol, comm, req + 4);
                }

                if (mype != 0)
                {
                    MPI_Irecv(&rbuf1[0], (m + 1) * m_loc, MPI_DOUBLE, destleft,
                        destleft, comm, req + 1);
                    MPI_Irecv(&actual[column[(n_loc - 1)]], 1, MPI_INT,
                        destleft, destleft + npecol, comm, req + 5);
                }

                if (mype != 0)
                {
                    MPI_Isend(&sbuf2[0], (m + 1) * m_loc, MPI_DOUBLE, destleft,
                        mype, comm, req + 2);
                    MPI_Isend(&scoll, 1, MPI_INT, destleft, mype + npecol, comm,
                        req + 6);
                    // cout<<"mype="<<mype<<", send col left
                    //"<<scoll<<endl;
                }

                if (mype != npecol - 1)
                {
                    MPI_Isend(&sbuf1[0], (m + 1) * m_loc, MPI_DOUBLE, destright,
                        mype, comm, req + 3);
                    MPI_Isend(&scolr, 1, MPI_INT, destright, mype + npecol,
                        comm, req + 7);
                    // cout<<"mype="<<mype<<", send col right
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
                        memcpy(&r_loc[k][column[(n_loc_2 - 1)] * m_loc],
                            &rbuf2[k * m_loc], m_loc * sizeof(double));
                    if (mype != 0)
                        memcpy(&r_loc[k][column[(n_loc - 1)] * m_loc],
                            &rbuf1[k * m_loc], m_loc * sizeof(double));
                }
                if (mype != npecol - 1)
                    memcpy(&u_loc[column[(n_loc_2 - 1)] * m_loc],
                        &rbuf2[m * m_loc], m_loc * sizeof(double));
                if (mype != 0)
                    memcpy(&u_loc[column[(n_loc - 1)] * m_loc],
                        &rbuf1[m * m_loc], m_loc * sizeof(double));
#endif

                // rotate pointers for local columns
                std::vector<int>::iterator bp = column.begin();
                if (mype == 0) bp++;
                std::vector<int>::iterator ep = column.end();
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
                    if (rp < m_loc)
                    {
                        sum += r_loc[k][rp + m_loc * ip]
                               * r_loc[k][rp + m_loc * ip];
                        for (int i = 0; i < m_loc; i++)
                            sum_off += r_loc[k][i + m_loc * ip]
                                       * r_loc[k][i + m_loc * ip];
                        sum_off -= r_loc[k][rp + m_loc * ip]
                                   * r_loc[k][rp + m_loc * ip];
                    }
                }
            }

#ifdef MGMOL_USE_SCALAPACK
            double in = sum;
            MPI_Allreduce(&in, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
            in = sum_off;
            MPI_Allreduce(&in, &sum_off, 1, MPI_DOUBLE, MPI_SUM, comm);
#endif
#if DEBUG
            std::cout << " Sum_off(A(k)):  " << sum_off << std::endl;
            std::cout << " Sum_diag(A(k)): " << sum << std::endl;
            std::cout << " Frobenius:        " << sum_off + sum << std::endl;
            std::cout << std::endl;
#endif

            isweep++;

            // delta = sum - oldsum;
            delta = std::abs(sum - oldsum);
#if DEBUG
            if (print_flag) std::cout << " sum:   " << sum << std::endl;
            if (print_flag) std::cout << " delta: " << delta << std::endl;
#endif

        } // isweep

        // loop over local columns
        for (int ip = 0; ip < n_loc; ip++)
        {
#if DEBUG
            // check columns back onto original PE
            if (!(actual[ip] < (mype + 1) * n_loc || actual[ip] == -1))
            {
                std::cout << "mype=" << mype << ", ip=" << ip
                          << ", actual[ip]=" << actual[ip] << std::endl;
            }
            if (!(actual[ip] >= mype * n_loc || actual[ip] == -1))
            {
                std::cout << "mype=" << mype << ", ip=" << ip
                          << ", actual[ip]=" << actual[ip] << std::endl;
            }
#endif
            assert(actual[ip] < (mype + 1) * n_loc || actual[ip] == -1);
            assert(actual[ip] >= mype * n_loc || actual[ip] == -1);

            if (actual[ip] < m_loc)
            {
                // update columns ip of r_
                for (int k = 0; k < m; k++)
                {
                    assert(ip < rmat[k]->nb());
                    rmat[k]->assignColumn(
                        &r_loc[k][m_loc * ip], actual[ip] - offset_);
                }
            }
        }

    } // if mype<npecol

#if DEBUG
    for (int k = 0; k < m; k++)
    {
        if (print_flag) std::cout << "rmat for k=" << k << std::endl;
        rmat[k]->print(cout, 0, 0, 5, 5);
    }
#endif

    // bcast u_loc into tmat
    for (int pe = 0; pe < npecol; pe++)
    {
        // bcast actual
        int root = pe;
        // loop over local columns
        for (int ip = 0; ip < n_loc; ip++)
        {
            int lip = actual[ip];
#ifdef MGMOL_USE_SCALAPACK
            mmpi.bcast(&lip, 1, root);
#endif
            if (lip < m_loc)
            {
                if (mype == root)
                {
                    memcpy(&tmat[m_loc * lip], &u_loc[m_loc * ip],
                        m_loc * sizeof(double));
                }
#ifdef MGMOL_USE_SCALAPACK
                mmpi.bcast(&tmat[m_loc * lip], m_loc, root);
#endif
            }
        }
    }

    if (print_flag)
    {
        std::cout << std::setprecision(12);
        std::cout << isweep << " sweeps in jade parallel" << std::endl;
        std::cout << " sum:   " << sum << std::endl;
        std::cout << " delta: " << delta << std::endl;
    }

    return delta;
}

#endif
