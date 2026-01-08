// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "GramMatrix.h"
#include "Power.h"
#include "ReplicatedMatrix.h"
#include "ReplicatedVector.h"

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrix.h"
#include "DistMatrix2SquareLocalMatrices.h"
#include "DistMatrixTools.h"
#include "DistVector.h"
#else
typedef double DISTMATDTYPE;
#endif

#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <vector>

template <class MatrixType>
GramMatrix<MatrixType>::GramMatrix(const int ndim) : dim_(ndim)
{
    matS_ = new MatrixType("S", ndim, ndim);
    ls_   = new MatrixType("LS", ndim, ndim);
    invS_ = new MatrixType("invS", ndim, ndim);

    work_ = new MatrixType("work", ndim, ndim);

    orbitals_index_ = -1;
    isLSuptodate_   = false;
    isInvSuptodate_ = false;
}

template <class MatrixType>
GramMatrix<MatrixType>::GramMatrix(const GramMatrix<MatrixType>& gm)
    : dim_(gm.dim_)
{
    matS_ = new MatrixType(*gm.matS_);
    ls_   = new MatrixType(*gm.ls_);
    invS_ = new MatrixType(*gm.invS_);

    work_ = new MatrixType(*gm.work_);

    orbitals_index_ = gm.orbitals_index_;
    isLSuptodate_   = gm.isLSuptodate_;
    isInvSuptodate_ = gm.isInvSuptodate_;
}

template <class MatrixType>
GramMatrix<MatrixType>& GramMatrix<MatrixType>::operator=(
    const GramMatrix<MatrixType>& gm)
{
    if (this == &gm) return *this;

    dim_ = gm.dim_;

    *matS_ = *gm.matS_;
    *ls_   = *gm.ls_;
    *invS_ = *gm.invS_;
    *work_ = *gm.work_;

    orbitals_index_ = gm.orbitals_index_;
    isLSuptodate_   = gm.isLSuptodate_;
    isInvSuptodate_ = gm.isInvSuptodate_;

    return *this;
}

template <class MatrixType>
GramMatrix<MatrixType>::~GramMatrix()
{
    delete matS_;
    delete invS_;
    delete ls_;
    delete work_;
}

template <class MatrixType>
void GramMatrix<MatrixType>::computeInverse()
{
    assert(isLSuptodate_);

    MatrixType zz(*ls_);
    zz.potri('l');
    work_->identity();
    invS_->symm('l', 'l', 1., zz, *work_, 0.);

    isInvSuptodate_ = true;
}

template <class MatrixType>
void GramMatrix<MatrixType>::transformLTML(
    MatrixType& mat, const DISTMATDTYPE alpha) const
{
    // mat = alpha* L**T * mat
    mat.trmm('l', 'l', 't', 'n', alpha, *ls_);
    // mat=mat*L
    mat.trmm('r', 'l', 'n', 'n', 1., *ls_);
}

// Solve Z<-L**(-T)*Z
template <class MatrixType>
void GramMatrix<MatrixType>::solveLST(MatrixType& z) const
{
    ls_->trtrs('l', 't', 'n', z);
}

#ifdef MGMOL_USE_SCALAPACK
template <>
double GramMatrix<dist_matrix::DistMatrix<double>>::computeCond()
{
#if 0
    double anorm   = matS_->norm('1');
    double invcond = ls_->pocon('l', anorm);
    double cond    = 0.;
    if (onpe0) cond = 1. / invcond;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.bcast(&cond, 1);
#else
    double emin;
    double emax;

    // static object so that eigenvectors inside Power
    // are saved from one call to the next
    static Power<dist_matrix::DistVector<double>,
        dist_matrix::DistMatrix<double>>
        power(dim_);

    power.computeEigenInterval(*matS_, emin, emax, 1.e-3);
    const double cond = emax / emin;
#endif

    return cond;
}
#endif

template <>
double GramMatrix<ReplicatedMatrix>::computeCond()
{
    const double cond = 1;

    return cond;
}

// mat is overwritten by inv(ls)*mat*inv(ls**T)
template <class MatrixType>
void GramMatrix<MatrixType>::sygst(MatrixType& mat)
{
    mat.sygst(1, 'l', *ls_);
}

template <class MatrixType>
void GramMatrix<MatrixType>::rotateAll(const MatrixType& matU)
{
    rotateSym(*matS_, matU, *work_);

    updateLS();

    computeInverse();
}

template <class MatrixType>
void GramMatrix<MatrixType>::setMatrix(
    const MatrixType& mat, const int orbitals_index)
{
    assert(matS_ != nullptr);

    *matS_          = mat;
    orbitals_index_ = orbitals_index;

    isInvSuptodate_ = false;

    updateLS();
}

template <class MatrixType>
void GramMatrix<MatrixType>::updateLS()
{
    // Cholesky decomposition of s
    *ls_     = *matS_;
    int info = ls_->potrf('l');
    if (info != 0)
    {
        std::cerr << "ERROR in GramMatrix<MatrixType>::updateLS()" << std::endl;
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.abort();
    }
    isLSuptodate_ = true;
}

template <class MatrixType>
void GramMatrix<MatrixType>::set2Id(const int orbitals_index)
{
    matS_->identity();
    invS_->identity();
    ls_->identity();

    orbitals_index_ = orbitals_index;
    isLSuptodate_   = true;
    isInvSuptodate_ = true;
}

template <class MatrixType>
double GramMatrix<MatrixType>::getLinDependent2states(int& st1, int& st2) const
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    std::vector<DISTMATDTYPE> eigenvalues(dim_);
    MatrixType u("u", dim_, dim_);
    // solve a standard symmetric eigenvalue problem
    MatrixType mat(*matS_);
    mat.syev('v', 'l', eigenvalues, u);
    // get the index of the two largest components of column 0 of u
    DISTMATDTYPE val1;
    st1 = u.iamax(0, val1);
    u.setVal(st1, 0, 0.);
    if (mmpi.instancePE0())
        std::cout
            << "GramMatrix<MatrixType>::getLinDependent2states... element val="
            << val1 << std::endl;
    DISTMATDTYPE val2;
    st2 = u.iamax(0, val2);
    if (mmpi.instancePE0())
        std::cout
            << "GramMatrix<MatrixType>::getLinDependent2states... element val="
            << val2 << std::endl;

    // look for 2nd largest coefficient with different sign
    while (val1 * val2 > 0.)
    {
        u.setVal(st2, 0, 0.);
        st2 = u.iamax(0, val2);
        if (mmpi.instancePE0())
            std::cout << "GramMatrix<MatrixType>::getLinDependent2states... "
                         "element val="
                      << val2 << std::endl;
    }

    return (double)eigenvalues[0];
}

template <class MatrixType>
double GramMatrix<MatrixType>::getLinDependent2states(
    int& st1, int& st2, int& st3) const
{
    std::vector<DISTMATDTYPE> eigenvalues(dim_);
    MatrixType u("u", dim_, dim_);
    // solve a standard symmetric eigenvalue problem
    MatrixType mat(*matS_);
    mat.syev('v', 'l', eigenvalues, u);
    // get the index of the two largest components of column 0 of u
    DISTMATDTYPE val1;
    st1 = u.iamax(0, val1);
    u.setVal(st1, 0, 0.);

    st2 = u.iamax(0, val1);
    u.setVal(st2, 0, 0.);

    st3 = u.iamax(0, val1);

    return (double)eigenvalues[0];
}

template <class MatrixType>
void GramMatrix<MatrixType>::computeLoewdinTransform(MatrixType& loewdinMat,
    std::shared_ptr<MatrixType> invLoewdin, const int orb_index)
{
    MatrixType mat(*matS_);
    MatrixType vect("vect", dim_, dim_);
    std::vector<DISTMATDTYPE> eigenvalues(dim_);
    mat.syev('v', 'l', eigenvalues, vect);

    std::vector<DISTMATDTYPE> diag_values(dim_);
    for (unsigned int i = 0; i < dim_; i++)
        diag_values[i] = (DISTMATDTYPE)(1. / sqrt(eigenvalues[i]));

    loewdinMat.clear();
    loewdinMat.setDiagonal(diag_values);

    //  loewdinMat = vect * D * vect^T, using mat as temporary storage
    mat.symm('r', 'l', 1., loewdinMat, vect, 0.);
    loewdinMat.gemm('n', 't', 1., mat, vect, 0.);

    if (invLoewdin) // compute inverse Loewdin matrix
    {
        // new Gram matrix is Identity
        set2Id(orb_index);

        for (unsigned int i = 0; i < dim_; i++)
            diag_values[i] = std::sqrt(eigenvalues[i]);
        invLoewdin->clear();
        invLoewdin->setDiagonal(diag_values);
        // invLoewdin = vect * D^-1 * vect^T, using mat as temporary storage
        mat.symm('r', 'l', 1., *invLoewdin, vect, 0.);
        invLoewdin->gemm('n', 't', 1., mat, vect, 0.);
    }
}

template <class MatrixType>
double GramMatrix<MatrixType>::getTraceDiagProductWithInvS(
    std::vector<DISTMATDTYPE>& ddiag)
{
    MatrixType diag("diag", dim_, dim_);
    diag.setDiagonal(ddiag);
    work_->gemm('n', 'n', 1., diag, *invS_, 0.);
    return work_->trace();
}

template <class MatrixType>
void GramMatrix<MatrixType>::applyInv(MatrixType& mat)
{
    assert(isLSuptodate_);

    ls_->potrs('l', mat);
}

template <class MatrixType>
template <class VectorType>
void GramMatrix<MatrixType>::applyInv(VectorType& mat)
{
    assert(isLSuptodate_);

    ls_->potrs('l', mat);
}

#ifdef MGMOL_USE_SCALAPACK
template class GramMatrix<dist_matrix::DistMatrix<DISTMATDTYPE>>;
template void GramMatrix<dist_matrix::DistMatrix<DISTMATDTYPE>>::applyInv(
    dist_matrix::DistVector<DISTMATDTYPE>&);
#endif
template class GramMatrix<ReplicatedMatrix>;
template void GramMatrix<ReplicatedMatrix>::applyInv(ReplicatedVector&);
