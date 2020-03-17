// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "GramMatrix.h"
#include "Control.h"
#include "DistMatrix2SquareLocalMatrices.h"
#include "DistVector.h"
#include "MGmol_MPI.h"
#include "Power.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <unistd.h>
#include <vector>

void rotateSym(dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
    const dist_matrix::DistMatrix<DISTMATDTYPE>& rotation_matrix,
    dist_matrix::DistMatrix<DISTMATDTYPE>& work);

GramMatrix::GramMatrix(const int ndim) : dim_(ndim)
{
    matS_ = new dist_matrix::DistMatrix<DISTMATDTYPE>("S", ndim, ndim);
    ls_   = new dist_matrix::DistMatrix<DISTMATDTYPE>("LS", ndim, ndim);
    invS_ = new dist_matrix::DistMatrix<DISTMATDTYPE>("invS", ndim, ndim);

    work_ = new dist_matrix::DistMatrix<DISTMATDTYPE>("work", ndim, ndim);

    orbitals_index_ = -1;
    isLSuptodate_   = false;
    isInvSuptodate_ = false;
}

GramMatrix::GramMatrix(const GramMatrix& gm) : dim_(gm.dim_)
{
    matS_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(*gm.matS_);
    ls_   = new dist_matrix::DistMatrix<DISTMATDTYPE>(*gm.ls_);
    invS_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(*gm.invS_);

    work_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(*gm.work_);

    orbitals_index_ = gm.orbitals_index_;
    isLSuptodate_   = gm.isLSuptodate_;
    isInvSuptodate_ = gm.isInvSuptodate_;
}

GramMatrix& GramMatrix::operator=(const GramMatrix& gm)
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

GramMatrix::~GramMatrix()
{
    delete matS_;
    delete invS_;
    delete ls_;
    delete work_;
}

void GramMatrix::computeInverse()
{
    assert(isLSuptodate_);

    dist_matrix::DistMatrix<DISTMATDTYPE> zz(*ls_);
    zz.potri('l');
    work_->identity();
    invS_->symm('l', 'l', 1., zz, *work_, 0.);

    isInvSuptodate_ = true;
}

void GramMatrix::transformLTML(
    dist_matrix::DistMatrix<DISTMATDTYPE>& mat, const DISTMATDTYPE alpha) const
{
    // mat = alpha* L**T * mat
    mat.trmm('l', 'l', 't', 'n', alpha, *ls_);
    // mat=mat*L
    mat.trmm('r', 'l', 'n', 'n', 1., *ls_);
}

// Solve Z<-L**(-T)*Z
void GramMatrix::solveLST(dist_matrix::DistMatrix<DISTMATDTYPE>& z) const
{
    ls_->trtrs('l', 't', 'n', z);
}

double GramMatrix::computeCond()
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

// mat is overwritten by inv(ls)*mat*inv(ls**T)
void GramMatrix::sygst(dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
{
    mat.sygst(1, 'l', *ls_);
}

void GramMatrix::rotateAll(const dist_matrix::DistMatrix<DISTMATDTYPE>& matU)
{
    rotateSym(*matS_, matU, *work_);

    updateLS();

    computeInverse();
}

void GramMatrix::setMatrix(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& mat, const int orbitals_index)
{
    assert(matS_ != nullptr);

    *matS_          = mat;
    orbitals_index_ = orbitals_index;

    isInvSuptodate_ = false;

    updateLS();
}

void GramMatrix::updateLS()
{
    // Cholesky decomposition of s
    *ls_     = *matS_;
    int info = ls_->potrf('l');
    if (info != 0)
    {
        print(*MPIdata::serr);
        if (onpe0)
            (*MPIdata::serr) << "ERROR in GramMatrix::updateLS()" << std::endl;
        sleep(5);
        Control& ct = *(Control::instance());
        ct.global_exit(2);
    }

    isLSuptodate_ = true;
}

void GramMatrix::set2Id(const int orbitals_index)
{
    matS_->identity();
    invS_->identity();
    ls_->identity();

    orbitals_index_ = orbitals_index;
    isLSuptodate_   = true;
    isInvSuptodate_ = true;
}

double GramMatrix::getLinDependent2states(int& st1, int& st2) const
{
    std::vector<DISTMATDTYPE> eigenvalues(dim_);
    dist_matrix::DistMatrix<DISTMATDTYPE> u("u", dim_, dim_);
    // solve a standard symmetric eigenvalue problem
    dist_matrix::DistMatrix<DISTMATDTYPE> mat(*matS_);
    mat.syev('v', 'l', eigenvalues, u);
    // get the index of the two largest components of column 0 of u
    DISTMATDTYPE val1;
    st1 = u.iamax(0, val1);
    u.setVal(st1, 0, 0.);
    if (onpe0)
        std::cout << "GramMatrix::getLinDependent2states... element val="
                  << val1 << std::endl;
    DISTMATDTYPE val2;
    st2 = u.iamax(0, val2);
    if (onpe0)
        std::cout << "GramMatrix::getLinDependent2states... element val="
                  << val2 << std::endl;

    // look for 2nd largest coefficient with different sign
    while (val1 * val2 > 0.)
    {
        u.setVal(st2, 0, 0.);
        st2 = u.iamax(0, val2);
        if (onpe0)
            std::cout << "GramMatrix::getLinDependent2states... element val="
                      << val2 << std::endl;
    }

    return (double)eigenvalues[0];
}

double GramMatrix::getLinDependent2states(int& st1, int& st2, int& st3) const
{
    std::vector<DISTMATDTYPE> eigenvalues(dim_);
    dist_matrix::DistMatrix<DISTMATDTYPE> u("u", dim_, dim_);
    // solve a standard symmetric eigenvalue problem
    dist_matrix::DistMatrix<DISTMATDTYPE> mat(*matS_);
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

void GramMatrix::computeLoewdinTransform(
    dist_matrix::DistMatrix<DISTMATDTYPE>& loewdinMat,
    std::shared_ptr<dist_matrix::DistMatrix<DISTMATDTYPE>> invLoewdin,
    std::shared_ptr<dist_matrix::DistMatrix<DISTMATDTYPE>> vect,
    const int orb_index)
{
    dist_matrix::DistMatrix<DISTMATDTYPE> mat(*matS_);
    std::vector<DISTMATDTYPE> eigenvalues(dim_);
    mat.syev('v', 'l', eigenvalues, *vect);

    std::vector<DISTMATDTYPE> diag_values(dim_);
    for (int i = 0; i < dim_; i++)
        diag_values[i] = (DISTMATDTYPE)(1. / sqrt(eigenvalues[i]));

    loewdinMat.clear();
    loewdinMat.setDiagonal(diag_values);

    //  loewdinMat = vect * D * vect^T, using mat as temporary storage
    mat.symm('r', 'l', 1., loewdinMat, *vect, 0.);
    loewdinMat.gemm('n', 't', 1., mat, *vect, 0.);

    if (invLoewdin) // compute inverse Loewdin matrix
    {
        // new Gram matrix is Identity
        set2Id(orb_index);

        for (unsigned int i = 0; i < dim_; i++)
            diag_values[i] = sqrt(eigenvalues[i]);
        invLoewdin->clear();
        invLoewdin->setDiagonal(diag_values);
        // invLoewdin = vect * D^-1 * vect^T, using mat as temporary storage
        mat.symm('r', 'l', 1., *invLoewdin, *vect, 0.);
        invLoewdin->gemm('n', 't', 1., mat, *vect, 0.);
    }
}
