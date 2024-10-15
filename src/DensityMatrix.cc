// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DensityMatrix.h"

#include "DistMatrix.h"
#include "MGmol_MPI.h"
#include "ReplicatedMatrix.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string.h>

const double factor_kernel4dot = 10.;

#define PROCRUSTES 0

// occupations in [0,1]
// DM eigenvalues in [0,orbital_occupation]

template <class MatrixType>
DensityMatrix<MatrixType>::DensityMatrix(const int ndim)
{
    assert(ndim > 0);

    dim_ = ndim;

    occ_uptodate_ = false;
    stripped_     = false;
    uniform_occ_  = false;

    MGmol_MPI& mmpi     = *(MGmol_MPI::instance());
    orbital_occupation_ = mmpi.nspin() > 1 ? 1. : 2.;

    orbitals_index_ = -1;

    dm_         = new MatrixType("DM", ndim, ndim);
    kernel4dot_ = new MatrixType("K4dot", ndim, ndim);
    work_       = new MatrixType("work", ndim, ndim);
    occupation_.resize(dim_);
    setDummyOcc();
}

template <class MatrixType>
DensityMatrix<MatrixType>::~DensityMatrix()
{
    assert(dm_ != nullptr);
    assert(kernel4dot_ != nullptr);
    assert(work_ != nullptr);

    delete dm_;
    delete kernel4dot_;
    delete work_;
}

template <class MatrixType>
void DensityMatrix<MatrixType>::build(const MatrixType& zmat,
    const std::vector<double>& occ, const int new_orbitals_index)
{
#ifdef PRINT_OPERATIONS
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (mmpi.instancePE0())
        std::cout << "template <class MatrixType> "
                     "DensityMatrix<MatrixType>::build(const "
                     "MatrixType&, const "
                     "vector<double>&, const int)"
                  << std::endl;
#endif

    // diagonal matrix with occ values in diagonal
    MatrixType gamma("Gamma", &occ[0], dim_, dim_);
    gamma.scal(orbital_occupation_); // rescale for spin

    // work_ = zmat*gamma with gamma symmetric
    work_->symm('r', 'l', 1., gamma, zmat, 0.);

    // dm_ = work_ * zmat^T
    dm_->gemm('n', 't', 1., *work_, zmat, 0.);

    std::vector<double> w(dim_);
    for (int i = 0; i < dim_; i++)
        w[i] = (double)(orbital_occupation_
                        * std::min(1., factor_kernel4dot * occ[i]));
    gamma.setDiagonal(w);

    work_->symm('r', 'l', 1., gamma, zmat, 0.);
    kernel4dot_->gemm('n', 't', 1., *work_, zmat, 0.);

    stripped_       = false;
    orbitals_index_ = new_orbitals_index;
}

template <class MatrixType>
void DensityMatrix<MatrixType>::build(
    const MatrixType& zmat, const int new_orbitals_index)
{
    build(zmat, occupation_, new_orbitals_index);
}

// build diagonal matrix
template <class MatrixType>
void DensityMatrix<MatrixType>::build(
    const std::vector<double>& occ, const int new_orbitals_index)
{
    assert(dm_ != nullptr);
    assert(!occ.empty());

    setOccupations(occ);

    build(new_orbitals_index);
}

template <class MatrixType>
void DensityMatrix<MatrixType>::build(const int new_orbitals_index)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
#ifdef PRINT_OPERATIONS
    if (mmpi.PE0())
        std::cout
            << "DensityMatrix<MatrixType>::build() for diagonal occupation..."
            << std::endl;
#endif
    if (!occ_uptodate_ && mmpi.instancePE0())
        std::cout << "Warning: occupations not up to date to build DM!!!"
                  << std::endl;

    MatrixType gamma("Gamma", &occupation_[0], dim_, dim_);
    gamma.scal(orbital_occupation_); // rescale for spin

    *dm_ = gamma;

    kernel4dot_->clear();
    std::vector<double> w(dim_);
    for (int i = 0; i < dim_; i++)
        w[i] = (double)(orbital_occupation_
                        * std::min(1., factor_kernel4dot * occupation_[i]));
    kernel4dot_->setDiagonal(w);

    stripped_       = false;
    orbitals_index_ = new_orbitals_index;
}

template <class MatrixType>
void DensityMatrix<MatrixType>::setUniform(
    const double nel, const int new_orbitals_index)
{
    assert(!occupation_.empty());

    MGmol_MPI& mmpi  = *(MGmol_MPI::instance());
    const double occ = (double)((double)nel / (double)dim_);
    if (mmpi.instancePE0())
        std::cout << "DensityMatrix::setUniform(), occupation = " << occ
                  << std::endl;
    assert(occ < 1.01);
    for (int i = 0; i < dim_; i++)
        occupation_[i] = occ;

    occ_uptodate_ = true;

    uniform_occ_ = true;

    build(occupation_, new_orbitals_index);
}

template <class MatrixType>
void DensityMatrix<MatrixType>::buildFromBlock(const MatrixType& block00)
{
    dm_->clear();
    dm_->assign(block00, 0, 0);
    dm_->print(std::cout, 0, 0, 25, 25);
}

template <class MatrixType>
void DensityMatrix<MatrixType>::rotate(
    const MatrixType& rotation_matrix, const bool flag_eigen)
{

    if (!flag_eigen)
    {
        MatrixType invU(rotation_matrix);
        std::vector<int> ipiv;
        invU.getrf(ipiv);

        // dm -> u**-1 * dm
        invU.getrs('n', *dm_, ipiv);

        // tmp = dm**T * u**-T
        MatrixType tmp(rotation_matrix);
        tmp.transpose(1., *dm_, 0.);

        // tmp = u**-1 * dm * u**-T
        invU.getrs('n', tmp, ipiv);

        *dm_ = tmp;
    }
}

template <class MatrixType>
void DensityMatrix<MatrixType>::printOccupations(std::ostream& os) const
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    if (mmpi.instancePE0())
    {
        os << std::endl << " Occupation numbers: ";

        // Print ten to a row.
        os.setf(std::ios::right, std::ios::adjustfield);
        os.setf(std::ios::fixed, std::ios::floatfield);
        os << std::setprecision(3);
        for (int i = 0; i < dim_; i++)
        {
            if ((i % 10) == 0) os << std::endl;
            os << std::setw(7) << occupation_[i] * orbital_occupation_ << " ";
        }

        os << std::endl;
    }
}

// double template <class MatrixType>
// DensityMatrix<MatrixType>::getSumOccupations()const
//{
//    double sum=0.;
//    for(int i=0;i<dim_;i++)
//    {
//        sum+=(double)occupation_[i]*(double)orbital_occupation_;
//    }
//
//    return sum;
//}

// solve the eigenvalue problem L^T*dm_*L*V=occ*V
// using the LL^T decomposition of S to get occ
template <class MatrixType>
void DensityMatrix<MatrixType>::diagonalize(
    const MatrixType& ls, std::vector<double>& occ)
{
    MatrixType evect("EigVect", dim_, dim_);

    *work_ = (*dm_);

    // *work_ = L**T * *work_
    work_->trmm('l', 'l', 't', 'n', 1., ls);
    // *work_ = *work_*L
    work_->trmm('r', 'l', 'n', 'n', 1., ls);

    // compute eigenvalues of work_
    work_->syev('n', 'l', occ, evect);
}

template <class MatrixType>
void DensityMatrix<MatrixType>::diagonalize(
    const char eigv, std::vector<double>& occ, MatrixType& vect)
{
    *work_ = (*dm_);

    // compute eigenvectors and eigenvalues of work_
    work_->syev(eigv, 'l', occ, vect);
}

template <class MatrixType>
void DensityMatrix<MatrixType>::computeOccupations(const MatrixType& ls)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
#ifdef PRINT_OPERATIONS
    if (mmpi.instancePE0())
        std::cout << "template <class MatrixType> "
                     "DensityMatrix<MatrixType>::computeOccupations()"
                  << std::endl;
#endif

    std::vector<double> occ(dim_);

    diagonalize(ls, occ);
    const double occinv = 1. / orbital_occupation_;

    const double tol = 1.e-5;
#ifndef NDEBUG
    const double tol_fail = 1.e-3;
#endif
    bool flag = false;
    for (int i = 0; i < dim_; i++)
    {
        double occ_val = (double)occ[i] * occinv;
        std::cout << std::setprecision(16);
        if (mmpi.instancePE0() && (occ_val < 0. - tol || occ_val > 1. + tol))
        {
            std::cout << "WARNING: template <class MatrixType> "
                         "DensityMatrix<MatrixType>::computeOccupations(), occ["
                      << i << "]=" << occ_val;
            // if( occ_uptodate_)std::cout<<" vs.
            // "<<occupation_[dim_-i-1];
            std::cout << std::endl;
            flag = true;
        }
        assert(occ_val > 0. - tol_fail);
        assert(occ_val < 1. + tol_fail);
        occ[i] = (double)std::max(0., occ_val);
        occ[i] = (double)std::min(1., occ_val);
    }
    if (flag) printOccupations(std::cout);

    for (int i = 0; i < dim_; i++)
    {
        occupation_[i] = occ[dim_ - i - 1];
    }
    occ_uptodate_ = true;
}

template <class MatrixType>
void DensityMatrix<MatrixType>::setOccupations(const std::vector<double>& occ)
{
    assert(!occ.empty());
#ifdef PRINT_OPERATIONS
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (mmpi.instancePE0())
        std::cout << "template <class MatrixType> "
                     "DensityMatrix<MatrixType>::setOccupations()"
                  << std::endl;
#endif
    assert((int)occ.size() == dim_);
    memcpy(&occupation_[0], &occ[0], dim_ * sizeof(double));
    occ_uptodate_ = true;
}

template <class MatrixType>
double DensityMatrix<MatrixType>::computeEntropy() const
{
    double s                  = 0.;
    const double tol          = 1.e-15;
    const double tol_interval = 1.e-6;

    assert(occ_uptodate_);
    for (int st = 0; st < dim_; st++)
    {
        const double fi = (double)occupation_[st];
        if (fi > 1. + tol_interval)
            std::cout << std::setprecision(15) << std::scientific << "f[" << st
                      << "]=" << fi << std::endl;
        assert(fi >= 0. - tol_interval);
        assert(fi <= 1. + tol_interval);
        if (fi < tol)
        {
            s += (1. - fi) * std::log(1. - fi);
        }
        else if (fi > 1. - tol)
        {
            s += fi * std::log(fi);
        }
        else
        {
            s += fi * std::log(fi) + (1. - fi) * std::log(1. - fi);
        }
    }

    return (double)(-orbital_occupation_) * s; // in units of kbt
}

template <class MatrixType>
void DensityMatrix<MatrixType>::setto2InvS(
    const MatrixType& invS, const int orbitals_index)
{
    *dm_ = invS;
    dm_->scal(orbital_occupation_);

    if (!occ_uptodate_)
    {
        for (int st = 0; st < dim_; st++)
            occupation_[st] = 1.;
        occ_uptodate_ = true;
    }
    uniform_occ_    = false;
    orbitals_index_ = orbitals_index;
}

template <class MatrixType>
void DensityMatrix<MatrixType>::stripS(const MatrixType& ls)
{
    assert(!stripped_);

    dm_->trmm('l', 'l', 't', 'n', 1., ls);
    dm_->trmm('r', 'l', 'n', 'n', 1., ls);

    uniform_occ_  = false;
    occ_uptodate_ = false;
    stripped_     = true;
}

template <class MatrixType>
void DensityMatrix<MatrixType>::dressUpS(
    const MatrixType& ls, const int new_orbitals_index)
{
    assert(stripped_);

    ls.trtrs('l', 't', 'n', *dm_);
    work_->transpose(1., *dm_, 0.);
    *dm_ = *work_;
    ls.trtrs('l', 't', 'n', *dm_);

    orbitals_index_ = new_orbitals_index;
    occ_uptodate_   = false;
    uniform_occ_    = false;
    stripped_       = false;
}

// dm_ -> u*dm_*u^T
// note: may lead to a dm with eigenvalues slightly outside [0.,1.]
// for u far from identity
template <class MatrixType>
void DensityMatrix<MatrixType>::transform(const MatrixType& u)
{
    work_->gemm('n', 't', 1., *dm_, u, 0.);
    dm_->gemm('n', 'n', 1., u, *work_, 0.);
}

template <class MatrixType>
double DensityMatrix<MatrixType>::getExpectation(const MatrixType& A)
{
    work_->gemm('n', 'n', 1., A, *dm_, 0.);
    return work_->trace();
}

template <class MatrixType>
void DensityMatrix<MatrixType>::mix(
    const double mix, const MatrixType& matA, const int new_orbitals_index)
{
    dm_->scal(1. - mix);

    dm_->axpy(mix, matA);
    orbitals_index_ = new_orbitals_index;
}

template class DensityMatrix<dist_matrix::DistMatrix<double>>;
#ifdef HAVE_MAGMA
template class DensityMatrix<ReplicatedMatrix>;
#endif
