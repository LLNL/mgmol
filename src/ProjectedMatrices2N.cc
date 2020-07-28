// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ProjectedMatrices2N.h"

template <class MatrixType>
ProjectedMatrices2N<MatrixType>::ProjectedMatrices2N(
    const int ndim, const bool with_spin)
    : ProjectedMatrices<MatrixType>(ndim, with_spin)
{
    bdim_ = ndim / 2;

    work2N_ = new MatrixType("work2N", ndim, ndim);
}

template <class MatrixType>
ProjectedMatrices2N<MatrixType>::~ProjectedMatrices2N()
{
    delete work2N_;
}

template <class MatrixType>
void ProjectedMatrices2N<MatrixType>::assignBlocksH(
    MatrixType& h11, MatrixType& h12, MatrixType& h21, MatrixType& h22)
{
    ProjectedMatrices<MatrixType>::matH_->assign(h11, 0, 0);
    ProjectedMatrices<MatrixType>::matH_->assign(h12, 0, bdim_);
    ProjectedMatrices<MatrixType>::matH_->assign(h21, bdim_, 0);
    ProjectedMatrices<MatrixType>::matH_->assign(h22, bdim_, bdim_);
}

template <class MatrixType>
void ProjectedMatrices2N<MatrixType>::iterativeUpdateDMwithEigenstates(
    const double occ_width, const int nel, const int iterative_index,
    const bool flag_reduce_T)
{
    const int dim = this->dim();
    std::vector<DISTMATDTYPE> eigenval(dim);

    ProjectedMatrices<MatrixType>::solveGenEigenProblem(*work2N_, eigenval);

    double kbT = occ_width;
    std::vector<DISTMATDTYPE> occ(dim);
    const DISTMATDTYPE tol = 1.e-6;
    double mu;
    do
    {
        if (onpe0)
            (*MPIdata::sout) << "MVP target with kbT = " << kbT << std::endl;
        mu = ProjectedMatrices<
            MatrixType>::computeChemicalPotentialAndOccupations(kbT, nel, dim);
        ProjectedMatrices<MatrixType>::getOccupations(occ);
        kbT *= 0.5;
    } while (occ[bdim_] > tol && flag_reduce_T);

    if (onpe0)
        (*MPIdata::sout) << "MVP target with mu = " << mu << " [Ry]"
                         << std::endl;
    ProjectedMatrices<MatrixType>::buildDM(*work2N_, iterative_index);
}

template class ProjectedMatrices2N<dist_matrix::DistMatrix<DISTMATDTYPE>>;
