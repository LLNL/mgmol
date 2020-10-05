// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cassert>
#include <vector>

#include "LinearSolver.h"
#include "LinearSolverMatrix.h"
#include "PreconILU.h"
#include "random.h"

Timer LinearSolver::gmres_tm_("LinearSolver::gmres");
Timer LinearSolver::eigmin_tm_("LinearSolver::computeEigmin");
Timer LinearSolver::eigmax_tm_("LinearSolver::computeEigmax");

LinearSolver::LinearSolver()
{
    // constructor
    iters_   = 0;
    resnorm_ = 0.0;
}
LinearSolver::~LinearSolver()
{
    // destructor
}

int LinearSolver::fgmres(const LinearSolverMatrix<lsdatatype>& LSMat,
    const PreconILU<pcdatatype>& precon, const int lrindex, double* sol,
    const double tol, const int im, const int maxits)
{
    int i, i1, ii, j, k, k1, its, im1, pti, pti1;
    int ptih = 0;
    int one  = 1;
    double *c, *s, *rs, t;
    double negt, eps1, gam, *vv, *z;
    double beta = 0.;

    const int n        = LSMat.n();
    const double coeff = -1.0;
    double scal        = 1.0;

    im1      = im + 1;
    int imsz = im1 * n;
    vv       = new double[(imsz)];
    memset(vv, 0, imsz * sizeof(double));
    imsz = im * n;
    z    = new double[(imsz)];
    memset(z, 0, imsz * sizeof(double));
    imsz = im1 * (im + 3);
    std::vector<double> hh(imsz, 0.);
    c  = hh.data() + im1 * im;
    s  = c + im1;
    rs = s + im1;
    /*-------------------- outer loop starts here */
    int retval = 0;
    its        = 0;
    memset(sol, 0, n * sizeof(double));
    if (LSMat.isrescaled()) scal = LSMat.getScale(lrindex);
    /*-------------------- Outer loop */
    while (its < maxits)
    {
        /*-------------------- compute initial residual vector */
        LSMat.matvec(sol, vv);
        DSCAL(&n, &coeff, vv, &one);

        vv[lrindex] = 1.0 * scal + vv[lrindex];
        beta        = DNRM2(&n, vv, &one);
        /*-------------------- print info --------- */
        //       if (onpe0 && fits != NULL && its == 0 && myid == 0 && DEBUG)
        //         printf("%8d   %10.2e\n",its, beta) ;

        if (beta == 0.0) break;
        t = 1.0 / beta;
        /*--------------------   normalize:  vv    =  vv   / beta */
        DSCAL(&n, &t, vv, &one);
        if (its == 0) eps1 = tol * beta;
        /*--------------------initialize 1-st term  of rhs of hessenberg mtx */
        rs[0] = beta;
        /*-------------------- Krylov loop*/
        i = -1;
        while ((i < im - 1) && (beta > eps1) && (its++ < maxits))
        {
            i++;
            i1   = i + 1;
            pti  = i * n;
            pti1 = i1 * n;
            /*------------------------------------------------------------
            |  (Right) Preconditioning Operation   z_{j} = M^{-1} v_{j}
            +-----------------------------------------------------------*/
            //          pc_solve_tm_.start();
            precon.LUsolve(&vv[pti], &z[pti]);
            // memcpy(&z[pti], &vv[pti],n*sizeof(double));
            // fgmres(LSMat, &vv[pti], &z[pti], 1.0e-3, im, maxits);
            //          pc_solve_tm_.stop();
            /*-------------------- matvec operation w = A z_{j} = A M^{-1} v_{j}
             */
            //          linear_solver_matvec_tm_.start();
            LSMat.matvec(&z[pti], &vv[pti1]);
            //          linear_solver_matvec_tm_.stop();
            /*-------------------- modified gram - schmidt...
            |     h_{i,j} = (w,v_{i});
            |     w  = w - h_{i,j} v_{i}
            +------------------------------------------------------------*/
            ptih = i * im1;
            for (j = 0; j <= i; j++)
            {
                t            = DDOT(&n, &vv[j * n], &one, &vv[pti1], &one);
                hh[ptih + j] = t;
                negt         = -t;
                DAXPY(&n, &negt, &vv[j * n], &one, &vv[pti1], &one);
            }
            /*-------------------- h_{j+1,j} = ||w||_{2}    */
            t             = DNRM2(&n, &vv[pti1], &one);
            hh[ptih + i1] = t;
            // report failure
            if (t == 0.0) return (1);
            t = 1.0 / t;
            /*-------------------- v_{j+1} = w / h_{j+1,j}  */
            DSCAL(&n, &t, &vv[pti1], &one);
            /*-------------------- done with modified gram schimdt/arnoldi step
            | now  update factorization of hh.
            | perform previous transformations  on i-th column of h
            +-------------------------------------------------------*/
            for (k = 1; k <= i; k++)
            {
                k1            = k - 1;
                t             = hh[ptih + k1];
                hh[ptih + k1] = c[k1] * t + s[k1] * hh[ptih + k];
                hh[ptih + k]  = -s[k1] * t + c[k1] * hh[ptih + k];
            }
            gam = std::sqrt(
                std::pow(hh[ptih + i], 2) + std::pow(hh[ptih + i1], 2));
            /*-------------------- check if gamma is zero */
            if (gam == 0.0) gam = epsmac;
            /*-------------------- get  next plane rotation    */
            c[i]   = hh[ptih + i] / gam;
            s[i]   = hh[ptih + i1] / gam;
            rs[i1] = -s[i] * rs[i];
            rs[i]  = c[i] * rs[i];
            /*-------------------- get residual norm + test convergence*/
            hh[ptih + i] = c[i] * hh[ptih + i] + s[i] * hh[ptih + i1];
            beta         = std::abs(rs[i1]);
            //              if(onpe0 && fits != NULL && myid == 0 && DEBUG)
            //	               fprintf(fits,"%8d   %10.2e\n", its, beta) ;
            /*-------------------- end second while loop [Arnoldi] */
        }
        /*-------------------- now compute solution. 1st, solve upper
                               triangular system*/
        rs[i] = rs[i] / hh[ptih + i];
        for (ii = i - 1; ii >= 0; ii--)
        {
            t = rs[ii];
            for (j = ii + 1; j <= i; j++)
                t -= hh[j * im1 + ii] * rs[j];
            rs[ii] = t / hh[ii * im1 + ii];
        }
        /*--------------------  linear combination of z_j's to get sol. */
        for (j = 0; j <= i; j++)
            DAXPY(&n, &rs[j], &z[j * n], &one, sol, &one);
        /*--------------------  restart outer loop if needed */
        if (beta < eps1)
            break;
        else
        {
            if (its >= maxits) retval = 1;
        }
        /*-------------------- end main while loop */
    }
    /*-------------------- prepare to return */
    iters_   = its;
    resnorm_ = beta;
    delete[] vv;
    delete[] z;
    return (retval);
}
/*-----------------end of fgmres --------------------------------------- */

int LinearSolver::fgmres(const LinearSolverMatrix<lsdatatype>& LSMat,
    const PreconILU<pcdatatype>& precon, const double* rhs, double* sol,
    const double tol, const int im, const int maxits)
{
    int i, i1, ii, j, k, k1, pti, pti1, ptih = 0;
    int one = 1;
    double t;
    double negt;
    double eps1 = 0.;
    double gam;
    double beta = 0.;

    const int n = LSMat.n();

    int im1  = im + 1;
    int imsz = im1 * n;
    std::unique_ptr<double[]> vv(new double[imsz]);
    memset(vv.get(), 0, imsz * sizeof(double));
    imsz = im * n;
    std::unique_ptr<double[]> z(new double[imsz]);
    memset(z.get(), 0, imsz * sizeof(double));
    imsz = im1 * (im + 3);
    std::unique_ptr<double[]> hh(new double[imsz]);
    memset(hh.get(), 0, imsz * sizeof(double));
    double* c  = hh.get() + im1 * im;
    double* s  = c + im1;
    double* rs = s + im1;
    /*-------------------- outer loop starts here */
    int retval = 0;
    int its    = 0;
    memset(sol, 0, n * sizeof(double));
    /*-------------------- Outer loop */
    while (its < maxits)
    {
        /*-------------------- compute initial residual vector */
        LSMat.matvec(sol, vv.get());
        for (j = 0; j < n; j++)
            vv[j] = rhs[j] - vv[j]; /*  vv[0]= initial residual */
        beta = DNRM2(&n, vv.get(), &one);
        /*-------------------- print info --------- */
        //       if (onpe0 && fits != NULL && its == 0 && myid == 0 && DEBUG)
        //         printf("%8d   %10.2e\n",its, beta) ;

        if (beta == 0.0) break;
        t = 1.0 / beta;
        /*--------------------   normalize:  vv    =  vv   / beta */
        DSCAL(&n, &t, vv.get(), &one);
        if (its == 0) eps1 = tol * beta;
        /*--------------------initialize 1-st term  of rhs of hessenberg mtx */
        rs[0] = beta;
        /*-------------------- Krylov loop*/
        i = -1;

        while ((i < im - 1) && (beta > eps1) && (its++ < maxits))
        {
            i++;
            i1   = i + 1;
            pti  = i * n;
            pti1 = i1 * n;
            /*------------------------------------------------------------
            |  (Right) Preconditioning Operation   z_{j} = M^{-1} v_{j}
            +-----------------------------------------------------------*/
            //          pc_solve_tm_.start();
            precon.LUsolve(&vv[pti], &z[pti]);
            // memcpy(&z[pti], &vv[pti],n*sizeof(double));
            // fgmres(LSMat, &vv[pti], &z[pti], 1.0e-3, im, maxits);
            //          pc_solve_tm_.stop();
            /*-------------------- matvec operation w = A z_{j} = A M^{-1} v_{j}
             */
            //          linear_solver_matvec_tm_.start();
            LSMat.matvec(&z[pti], &vv[pti1]);
            //          linear_solver_matvec_tm_.stop();
            /*-------------------- modified gram - schmidt...
            |     h_{i,j} = (w,v_{i});
            |     w  = w - h_{i,j} v_{i}
            +------------------------------------------------------------*/
            ptih = i * im1;
            for (j = 0; j <= i; j++)
            {
                t            = DDOT(&n, &vv[j * n], &one, &vv[pti1], &one);
                hh[ptih + j] = t;
                negt         = -t;
                DAXPY(&n, &negt, &vv[j * n], &one, &vv[pti1], &one);
            }

            /*-------------------- h_{j+1,j} = ||w||_{2}    */
            t             = DNRM2(&n, &vv[pti1], &one);
            hh[ptih + i1] = t;

            // report failure
            if (t == 0.0) return (1);
            t = 1.0 / t;
            /*-------------------- v_{j+1} = w / h_{j+1,j}  */
            DSCAL(&n, &t, &vv[pti1], &one);
            /*-------------------- done with modified gram schimdt/arnoldi step
            | now  update factorization of hh.
            | perform previous transformations  on i-th column of h
            +-------------------------------------------------------*/
            for (k = 1; k <= i; k++)
            {
                k1            = k - 1;
                t             = hh[ptih + k1];
                hh[ptih + k1] = c[k1] * t + s[k1] * hh[ptih + k];
                hh[ptih + k]  = -s[k1] * t + c[k1] * hh[ptih + k];
            }

            gam = std::sqrt(
                std::pow(hh[ptih + i], 2) + std::pow(hh[ptih + i1], 2));
            /*-------------------- check if gamma is zero */
            if (gam == 0.0) gam = epsmac;
            /*-------------------- get  next plane rotation    */
            c[i]   = hh[ptih + i] / gam;
            s[i]   = hh[ptih + i1] / gam;
            rs[i1] = -s[i] * rs[i];
            rs[i]  = c[i] * rs[i];
            /*-------------------- get residual norm + test convergence*/
            hh[ptih + i] = c[i] * hh[ptih + i] + s[i] * hh[ptih + i1];
            beta         = std::abs(rs[i1]);
            //              if(onpe0 && fits != NULL && myid == 0 && DEBUG)
            //	               fprintf(fits,"%8d   %10.2e\n", its, beta) ;
            /*-------------------- end second while loop [Arnoldi] */
        }
        /*-------------------- now compute solution. 1st, solve upper
                               triangular system*/
        rs[i] = rs[i] / hh[ptih + i];
        for (ii = i - 1; ii >= 0; ii--)
        {
            t = rs[ii];
            for (j = ii + 1; j <= i; j++)
                t -= hh[j * im1 + ii] * rs[j];
            rs[ii] = t / hh[ii * im1 + ii];
        }
        /*--------------------  linear combination of z_j's to get sol. */
        for (j = 0; j <= i; j++)
            DAXPY(&n, &rs[j], &z[j * n], &one, sol, &one);
        /*--------------------  restart outer loop if needed */
        if (beta < eps1)
            break;
        else
        {
            if (its >= maxits) retval = 1;
        }
        /*-------------------- end main while loop */
    }

    /*-------------------- prepare to return */
    iters_   = its;
    resnorm_ = beta;
    return (retval);
}
/*-----------------end of fgmres ---------------------------------------*/

// Compute minimum eigenvalue of LSMat - Assumes SPD matrix (requires Rayleigh
// quotient to retrieve sign of eigenvalue)
double LinearSolver::computeEigMin(
    LinearSolverMatrix<lsdatatype>& LSMat, std::vector<double>& evec)
{
    int one = 1, n = LSMat.n(), maxits = 10;
    double gamma = 0., beta = 0., eigmin;

    eigmin_tm_.start();

    double* rhs = &evec[0];

    // initialize sol vector
    double* solptr = new double[n];
    memset(solptr, 0, n * sizeof(double));

    // compute preconditioner -- ILU0
    PreconILU<pcdatatype> precon(LSMat, 1.0e-3, 1000, 0);
    /* setup the preconditioner */
    precon.setup(LSMat, PCILU0);
    // begin
    if (n == 1) maxits = 1; // just do one iteration for trivial case
    for (int i = 0; i < maxits; i++)
    {
        // solve
        fgmres(LSMat, precon, rhs, solptr, 1.0e-6, 10, 10);

        gamma  = DNRM2(&n, solptr, &one);
        eigmin = 1. / gamma;
        // normalize solution
        DSCAL(&n, &eigmin, solptr, &one);
        // update evec -- set rhs = solptr
        memcpy(rhs, solptr, n * sizeof(double));
        // check for convergence
        beta -= eigmin;
        if (std::abs(beta) < 1.0e-3) break;

        // reset beta
        beta = eigmin;
    }

    delete[] solptr;

    eigmin_tm_.stop();
    return (eigmin);
}

// Compute maximum eigenvalue of LSMat - Assumes SPD matrix (requires Rayleigh
// quotient to retrieve sign of eigenvalue)
double LinearSolver::computeEigMax(
    LinearSolverMatrix<lsdatatype>& LSMat, std::vector<double>& evec)
{
    int one = 1, n = LSMat.n(), maxits = 10;
    double gamma = 0., beta = 0., eigmax;

    eigmax_tm_.start();

    // get pointer to evec
    std::vector<double>::iterator it = evec.begin();
    double* rhs                      = &(*it);
    // initialize sol vector
    double* solptr = new double[n];
    memset(solptr, 0, n * sizeof(double));

    for (int i = 0; i < maxits; i++)
    {
        // do matvec
        LSMat.matvec(rhs, solptr);
        // compute norm
        eigmax = DNRM2(&n, solptr, &one);
        // normalize result
        gamma = 1. / eigmax;
        DSCAL(&n, &gamma, solptr, &one);
        // update evec
        memcpy(rhs, solptr, n * sizeof(double));
        // check for convergence
        beta -= eigmax;
        if (std::abs(beta) < 1.0e-3) break;

        // reset beta
        beta = eigmax;
    }

    delete[] solptr;

    eigmax_tm_.stop();
    return (eigmax);
}

// Compute condition number of LSMat
double LinearSolver::computeConditionNumber(
    LinearSolverMatrix<lsdatatype>& LSMat)
{
    assert(LSMat.n() > 0);

    srand(13579);

    const int n    = LSMat.n();
    double condest = 0.;
    //   const double nrm = 1./(double)n;

    //   std::vector<double> randvec(generate_rand(n));
    std::vector<double> minvec(generate_rand(n));
    std::vector<double> maxvec(minvec);

    double emin = computeEigMin(LSMat, minvec);
    double emax = computeEigMax(LSMat, maxvec);

    condest = emax / emin;

    return condest;
}
/* mixed precision fgmres */
/* vector norms and dot products are returned in double precision (used in the
 * MG-S/ Arnoldi process). vectors and matrices used to perform transformations
 * on the Hessenberg matrix are stored in single precision. Updates to the
 * solution are stored in single precision, but the updates themselves are done
 * in double precision. Matrix and vector operations (linear algebra ops) on
 * single precision data are 'computed' in double precision and the results are
 * stored in single precision
 */
int LinearSolver::fgmres_mp(LinearSolverMatrix<lsdatatype>& LSMat,
    PreconILU<pcdatatype>& precon, const int lrindex, double* sol,
    const double tol, const int im, const int maxits)
{
    int i, i1, ii, j, k, k1, its, im1, pti, pti1, ptih = 0;

    double eps1, gam, negt, t;
    double beta = 0.;

    const int n        = LSMat.n();
    const double coeff = -1.0;
    double scal        = 1.0;

    float *hh, *c, *s, *rs, *vv, *z;
    im1              = im + 1;
    int vvsz         = im1 * n;
    int zsz          = im * n;
    int hhsz         = im1 * im;
    int gsz          = 3 * im1;
    int storage_size = vvsz + zsz + hhsz + gsz;
    vv               = new float[storage_size];
    memset(vv, 0, storage_size * sizeof(float));
    z  = vv + vvsz;
    hh = z + zsz;
    c  = hh + hhsz;
    s  = c + im1;
    rs = s + im1;

    //    double *ptrv = new double[n];
    //    double *ptrz = new double[zsz];
    //    memset(ptrz, 0, zsz*sizeof(double));

    /*-------------------- outer loop starts here */
    int retval = 0;
    its        = 0;
    memset(sol, 0, n * sizeof(double));
    if (LSMat.isrescaled()) scal = LSMat.getScale(lrindex);
    /*-------------------- Outer loop */
    while (its < maxits)
    {
        /*-------------------- compute initial residual vector */
        LSMat.matvec(sol, vv);
        LinearAlgebraUtils<MemorySpace::Host>::MPscal(n, coeff, vv);

        vv[lrindex] = 1.0 * scal + vv[lrindex];
        beta        = MPnrm2(n, vv);
        /*-------------------- print info --------- */
        //       if (onpe0 && fits != NULL && its == 0 && myid == 0 && DEBUG)
        //         printf("%8d   %10.2e\n",its, beta) ;

        if (beta == 0.0) break;
        t = 1.0 / beta;
        /*--------------------   normalize:  vv    =  vv   / beta */
        LinearAlgebraUtils<MemorySpace::Host>::MPscal(n, t, vv);
        if (its == 0) eps1 = tol * beta;
        /*--------------------initialize 1-st term  of rhs of hessenberg mtx */
        rs[0] = (float)beta;
        /*-------------------- Krylov loop*/
        i = -1;
        while ((i < im - 1) && (beta > eps1) && (its++ < maxits))
        {
            i++;
            i1   = i + 1;
            pti  = i * n;
            pti1 = i1 * n;
            /*------------------------------------------------------------
            |  (Right) Preconditioning Operation   z_{j} = M^{-1} v_{j}
            +-----------------------------------------------------------*/
            /*
                      for(int ik=0; ik<n; ik++)
                      {
                         ptrv[ik] = (double)vv[pti+ik];
                         ptrz[pti+ik] = (double)z[pti+ik];
                      }
                      precon.LUsolve(ptrv, &ptrz[pti]);
                      for(int ik=0; ik<n; ik++)
                      {
                         vv[pti+ik] = (float)ptrv[ik];
                         z[pti+ik] = (float)ptrz[pti+ik];
                      }
            */
            precon.LUsolve(&vv[pti], &z[pti]);

            /*-------------------- matvec operation w = A z_{j} = A M^{-1} v_{j}
             */
            //          LSMat.matvec(&ptrz[pti], &vv[pti1]);
            LSMat.matvec(&z[pti], &vv[pti1]);
            /*-------------------- modified gram - schmidt...
            |     h_{i,j} = (w,v_{i});
            |     w  = w - h_{i,j} v_{i}
            +------------------------------------------------------------*/
            ptih = i * im1;
            for (j = 0; j <= i; j++)
            {
                t = LinearAlgebraUtils<MemorySpace::Host>::MPdot(
                    n, &vv[j * n], &vv[pti1]);
                hh[ptih + j] = static_cast<float>(t);
                negt         = -t;
                LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
                    n, negt, &vv[j * n], &vv[pti1]);
            }
            /*-------------------- h_{j+1,j} = ||w||_{2}    */
            t             = MPnrm2(n, &vv[pti1]);
            hh[ptih + i1] = static_cast<float>(t);
            // report failure
            if (t == 0.0) return (1);
            t = 1.0 / t;
            /*-------------------- v_{j+1} = w / h_{j+1,j}  */
            LinearAlgebraUtils<MemorySpace::Host>::MPscal(n, t, &vv[pti1]);
            /*-------------------- done with modified gram schimdt/arnoldi step
            | now  update factorization of hh.
            | perform previous transformations  on i-th column of h
            +-------------------------------------------------------*/
            for (k = 1; k <= i; k++)
            {
                k1            = k - 1;
                t             = (double)hh[ptih + k1];
                hh[ptih + k1] = (float)((double)c[k1] * t
                                        + (double)s[k1] * (double)hh[ptih + k]);
                hh[ptih + k]  = (float)((double)-s[k1] * t
                                       + (double)c[k1] * (double)hh[ptih + k]);
            }
            gam = std::sqrt(std::pow((double)hh[ptih + i], 2)
                            + std::pow((double)hh[ptih + i1], 2));
            /*-------------------- check if gamma is zero */
            if (gam == 0.0) gam = epsmac;
            /*-------------------- get  next plane rotation    */
            c[i]   = (float)((double)hh[ptih + i] / gam);
            s[i]   = (float)((double)hh[ptih + i1] / gam);
            rs[i1] = (float)((double)-s[i] * (double)rs[i]);
            rs[i]  = (float)((double)c[i] * (double)rs[i]);
            /*-------------------- get residual norm + test convergence*/
            hh[ptih + i] = (float)((double)c[i] * (double)hh[ptih + i]
                                   + (double)s[i] * (double)hh[ptih + i1]);
            beta         = std::abs((double)rs[i1]);
            //              if(onpe0 && fits != NULL && myid == 0 && DEBUG)
            //	               fprintf(fits,"%8d   %10.2e\n", its, beta) ;
            /*-------------------- end second while loop [Arnoldi] */
        }
        /*-------------------- now compute solution. 1st, solve upper
                               triangular system*/
        rs[i] = (float)((double)rs[i] / (double)hh[ptih + i]);
        for (ii = i - 1; ii >= 0; ii--)
        {
            t = (double)rs[ii];
            for (j = ii + 1; j <= i; j++)
                t -= (double)hh[j * im1 + ii] * (double)rs[j];
            rs[ii] = (float)(t / (double)hh[ii * im1 + ii]);
        }
        /*--------------------  linear combination of z_j's to get sol. */
        for (j = 0; j <= i; j++)
        {
            LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
                n, rs[j], &z[j * n], sol);
        }
        /*--------------------  restart outer loop if needed */
        if (beta < eps1)
            break;
        else
        {
            if (its >= maxits) retval = 1;
        }
        /*-------------------- end main while loop */
    }
    /*-------------------- prepare to return */
    iters_   = its;
    resnorm_ = beta;
    delete[] vv;
    //  delete [] ptrz;
    //  delete [] ptrv;
    return (retval);
}
/*-----------------end of fgmres --------------------------------------- */
