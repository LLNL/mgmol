// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ProjectedMatrices2N.h"

ProjectedMatrices2N::ProjectedMatrices2N(const int ndim, const bool with_spin)
    : ProjectedMatrices(ndim, with_spin)
{
    bdim_ = ndim / 2;

    work2N_ = new dist_matrix::DistMatrix<DISTMATDTYPE>("work2N", ndim, ndim);
}

ProjectedMatrices2N::~ProjectedMatrices2N() { delete work2N_; }

void ProjectedMatrices2N::assignBlocksH(
    dist_matrix::DistMatrix<DISTMATDTYPE>& h11,
    dist_matrix::DistMatrix<DISTMATDTYPE>& h12,
    dist_matrix::DistMatrix<DISTMATDTYPE>& h21,
    dist_matrix::DistMatrix<DISTMATDTYPE>& h22)
{
    matH_->assign(h11, 0, 0);
    matH_->assign(h12, 0, bdim_);
    matH_->assign(h21, bdim_, 0);
    matH_->assign(h22, bdim_, bdim_);
}

void ProjectedMatrices2N::iterativeUpdateDMwithEigenstates(
    const double occ_width, const int nel, const int iterative_index,
    const bool flag_reduce_T)
{
    const int dim = this->dim();
    vector<DISTMATDTYPE> eigenval(dim);

    solveGenEigenProblem(*work2N_, eigenval);
    setAuxilliaryEnergiesFromEigenenergies();

    double kbT = occ_width;
    vector<DISTMATDTYPE> occ(dim);
    const DISTMATDTYPE tol = 1.e-6;
    double mu;
    do
    {
        if (onpe0) (*MPIdata::sout) << "MVP target with kbT = " << kbT << endl;
        mu = computeChemicalPotentialAndOccupations(kbT, nel, dim);
        getOccupations(occ);
        kbT *= 0.5;
    } while (occ[bdim_] > tol && flag_reduce_T);

    if (onpe0)
        (*MPIdata::sout) << "MVP target with mu = " << mu << " [Ry]" << endl;
    buildDM(*work2N_, iterative_index);
}
