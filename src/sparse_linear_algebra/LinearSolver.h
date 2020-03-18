// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
 * Class to handle linear system operations
 */

#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_

#include "LinearSolverMatrix.h"
#include "PreconILU.h"
#include "mputils.h"
#include <vector>

#define epsmac 1.0e-16

class LinearSolver
{
    static Timer gmres_tm_;
    static Timer eigmin_tm_;
    static Timer eigmax_tm_;

    int n_; /* matrix size */
    short iters_; // fgmres iteration count
    double resnorm_; // fgmres residual norm

    // gmres no rhs vector
    int fgmres(const LinearSolverMatrix<lsdatatype>& LSMat,
        const PreconILU<pcdatatype>& precon, const int lrindex, double* sol,
        const double tol, const int im, const int maxits);
    // mixed precision variant
    int fgmres_mp(LinearSolverMatrix<lsdatatype>& LSMat,
        PreconILU<pcdatatype>& precon, const int lrindex, double* sol,
        const double tol, const int im, const int maxits);
    // gmres with rhs vector
    int fgmres(const LinearSolverMatrix<lsdatatype>& LSMat,
        const PreconILU<pcdatatype>& precon, const double* rhs, double* sol,
        const double tol, const int im, const int maxits);

public:
    LinearSolver(); // default constructor
    LinearSolver(LinearSolverMatrix<lsdatatype> LSMat);
    // compute minimum eigenvalue of LSMat
    double computeEigMin(
        LinearSolverMatrix<lsdatatype>& LSMat, std::vector<double>& evec);
    // compute maximum eigenvalue of LSMat
    double computeEigMax(
        LinearSolverMatrix<lsdatatype>& LSMat, std::vector<double>& evec);
    // compute condition number of LSMat
    double computeConditionNumber(LinearSolverMatrix<lsdatatype>& LSMat);
    ~LinearSolver(); // destructor

    // solve the linear system
    int solve(const LinearSolverMatrix<lsdatatype>& LSMat,
        const PreconILU<pcdatatype>& precon, const int lrindex, double* sol,
        const double tol, const int im, const int maxits)
    {
        gmres_tm_.start();
        int conv = fgmres(LSMat, precon, lrindex, sol, tol, im, maxits);
        //       int conv = fgmres_mp(LSMat, precon, lrindex, sol, tol, im,
        //       maxits);
        gmres_tm_.stop();
        return conv;
    }
    int solve(const LinearSolverMatrix<lsdatatype>& LSMat,
        const PreconILU<pcdatatype>& precon, const double* rhs, double* sol,
        const double tol, const int im, const int maxits)
    {
        gmres_tm_.start();
        int conv = fgmres(LSMat, precon, rhs, sol, tol, im, maxits);
        gmres_tm_.stop();
        return conv;
    }
    // get iteration count
    short iters() { return (iters_); }
    // get residual norm
    double residualNorm() { return (resnorm_); }
    // print timers
    static void printTimers(std::ostream& os)
    {
        gmres_tm_.print(os);
        eigmin_tm_.print(os);
    }
};

#endif
