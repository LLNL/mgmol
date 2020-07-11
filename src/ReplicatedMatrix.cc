// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifdef HAVE_MAGMA

#include "ReplicatedMatrix.h"
#include "memory_space.h"
#include "random.h"

#include "magma_v2.h"

#include <iostream>

using MemoryDev = MemorySpace::Memory<double, MemorySpace::Device>;

constexpr double gpuroundup = 32;

void rotateSym(ReplicatedMatrix& mat, const ReplicatedMatrix& rotation_matrix,
    ReplicatedMatrix& work)
{
    work.symm('l', 'l', 1., mat, rotation_matrix, 0.);
    mat.gemm('t', 'n', 1., rotation_matrix, work, 0.);
}

ReplicatedMatrix::ReplicatedMatrix(
    const std::string name, const int m, const int n)
    : dim_(m),
      ld_(magma_roundup(dim_, gpuroundup)),
      device_data_(MemoryDev::allocate(dim_ * ld_), MemoryDev::free)
{
    assert(m == n);
}

ReplicatedMatrix::ReplicatedMatrix(const std::string name, const int n)
    : dim_(n),
      ld_(magma_roundup(dim_, gpuroundup)),
      device_data_(MemoryDev::allocate(dim_ * ld_), MemoryDev::free)
{
}

ReplicatedMatrix::ReplicatedMatrix(const ReplicatedMatrix& mat)
    : dim_(mat.dim_),
      ld_(mat.ld_),
      device_data_(MemoryDev::allocate(dim_ * ld_), MemoryDev::free)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dcopymatrix(dim_, dim_, mat.device_data_.get(), mat.ld_,
        device_data_.get(), ld_, magma_singleton.queue_);
}

ReplicatedMatrix& ReplicatedMatrix::operator=(const ReplicatedMatrix& rhs)
{
    if (this != &rhs)
    {
        ld_  = rhs.ld_;
        dim_ = rhs.dim_;
        device_data_.reset(MemoryDev::allocate(dim_ * ld_));

        auto& magma_singleton = MagmaSingleton::get_magma_singleton();

        magma_dcopymatrix(dim_, dim_, rhs.device_data_.get(), rhs.ld_,
            device_data_.get(), ld_, magma_singleton.queue_);
    }
    return *this;
}

ReplicatedMatrix::~ReplicatedMatrix() {}

void ReplicatedMatrix::axpy(const double alpha, const ReplicatedMatrix& a)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magmablas_dgeadd(dim_, dim_, alpha, a.device_data_.get(), a.ld_,
        device_data_.get(), ld_, magma_singleton.queue_);
}

void ReplicatedMatrix::setRandom(const double minv, const double maxv)
{
    std::vector<double> mat(dim_ * dim_);

    generateRandomData(mat, minv, maxv);

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dsetmatrix(dim_, dim_, mat.data(), dim_, device_data_.get(), ld_,
        magma_singleton.queue_);
}

void ReplicatedMatrix::identity()
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magmablas_dlaset(MagmaFull, dim_, dim_, 0.0, 1.0, device_data_.get(), ld_,
        magma_singleton.queue_);
}

// this = alpha * transpose(A) + beta * this
void ReplicatedMatrix::transpose(
    const double alpha, const ReplicatedMatrix& a, const double beta)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    double* dwork;
    magma_int_t ret = magma_dmalloc(&dwork, dim_ * ld_);
    if (ret != MAGMA_SUCCESS)
    {
        std::cerr << "magma_dmalloc failed!" << std::endl;
    }

    magmablas_dtranspose(dim_, dim_, a.device_data_.get(), a.ld_, dwork, ld_,
        magma_singleton.queue_);

    magmablas_dgeadd2(dim_, dim_, alpha, dwork, ld_, beta, device_data_.get(),
        ld_, magma_singleton.queue_);

    magma_singleton.sync();
    magma_free(dwork);
}

void ReplicatedMatrix::gemm(const char transa, const char transb,
    const double alpha, const ReplicatedMatrix& a, const ReplicatedMatrix& b,
    const double beta)
{
    magma_trans_t magma_transa = magma_trans_const(transa);
    magma_trans_t magma_transb = magma_trans_const(transb);

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magmablas_dgemm(magma_transa, magma_transb, dim_, dim_, dim_, alpha,
        a.device_data_.get(), a.ld_, b.device_data_.get(), b.ld_, beta,
        device_data_.get(), ld_, magma_singleton.queue_);
}

void ReplicatedMatrix::symm(const char side, const char uplo,
    const double alpha, const ReplicatedMatrix& a, const ReplicatedMatrix& b,
    const double beta)
{
    magma_side_t magma_side = magma_side_const(side);
    magma_uplo_t magma_uplo = magma_uplo_const(uplo);

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dsymm(magma_side, magma_uplo, dim_, dim_, alpha, a.device_data_.get(),
        a.ld_, b.device_data_.get(), b.ld_, beta, device_data_.get(), ld_,
        magma_singleton.queue_);
}

int ReplicatedMatrix::potrf(char uplo)
{
    magma_uplo_t magma_uplo = magma_uplo_const(uplo);

    int info;
    magma_dpotrf_gpu(magma_uplo, dim_, device_data_.get(), ld_, &info);
    if (info != 0)
        std::cerr << "magma_dpotrf_gpu failed, info = " << info << std::endl;

    return info;
}

int ReplicatedMatrix::potri(char uplo)
{
    magma_uplo_t magma_uplo = magma_uplo_const(uplo);

    int info;
    magma_dpotri_gpu(magma_uplo, dim_, device_data_.get(), ld_, &info);
    if (info != 0)
        std::cerr << "magma_dpotri_gpu failed, info = " << info << std::endl;

    return info;
}

// Solve a system of linear equations A*X = B with a symmetric
// positive definite matrix A using the Cholesky factorization
// A = U**T*U or A = L*L**T computed by potrf
void ReplicatedMatrix::potrs(char uplo, ReplicatedMatrix& b)
{
    magma_uplo_t magma_uplo = magma_uplo_const(uplo);

    int info;
    magma_dpotrs_gpu(magma_uplo, dim_, dim_, device_data_.get(), ld_,
        b.device_data_.get(), b.ld_, &info);
    if (info != 0)
        std::cerr << "magma_dpotrs_gpu failed, info = " << info << std::endl;
}

void ReplicatedMatrix::syev(
    char jobz, char uplo, std::vector<double>& evals, ReplicatedMatrix& z)
{
    magma_vec_t magma_jobz  = magma_vec_const(jobz);
    magma_uplo_t magma_uplo = magma_uplo_const(uplo);

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    // copy matrix into z
    magmablas_dlacpy(MagmaFull, dim_, dim_, device_data_.get(), ld_,
        z.device_data_.get(), z.ld_, magma_singleton.queue_);
    magma_int_t nb = magma_get_ssytrd_nb(dim_);
    magma_int_t lwork
        = std::max(2 * dim_ + dim_ * nb, 1 + 6 * dim_ + 2 * dim_ * dim_);
    int liwork = 3 + 5 * dim_;

    int info;
    std::vector<double> wa(dim_ * dim_);
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);

    magma_dsyevd_gpu(magma_jobz, magma_uplo, dim_, z.device_data_.get(), z.ld_,
        evals.data(), wa.data(), dim_, work.data(), lwork, iwork.data(), liwork,
        &info);
    if (info != 0)
        std::cerr << "magma_dsyevd_gpu failed, info = " << info << std::endl;
    // for(auto& d : evals)std::cout<<d<<std::endl;
}

void ReplicatedMatrix::sygst(int itype, char uplo, const ReplicatedMatrix& b)
{
    magma_uplo_t magma_uplo = magma_uplo_const(uplo);
    magma_int_t magma_itype = static_cast<magma_int_t>(itype);

    int info;
    magma_dsygst_gpu(magma_itype, magma_uplo, dim_, device_data_.get(), ld_,
        b.device_data_.get(), ld_, &info);
    if (info != 0)
        std::cerr << "magma_dsygst_gpu failed, info = " << info << std::endl;
}

void ReplicatedMatrix::trmm(const char side, const char uplo, const char trans,
    const char diag, const double alpha, const ReplicatedMatrix& a)
{
    magma_side_t magma_side   = magma_side_const(side);
    magma_uplo_t magma_uplo   = magma_uplo_const(uplo);
    magma_trans_t magma_trans = magma_trans_const(trans);
    magma_diag_t magma_diag   = magma_diag_const(diag);

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dtrmm(magma_side, magma_uplo, magma_trans, magma_diag, dim_, dim_,
        alpha, a.device_data_.get(), a.ld_, device_data_.get(), ld_,
        magma_singleton.queue_);
}

void ReplicatedMatrix::trtrs(const char uplo, const char trans, const char diag,
    ReplicatedMatrix& b) const
{
    magma_uplo_t magma_uplo   = magma_uplo_const(uplo);
    magma_trans_t magma_trans = magma_trans_const(trans);
    magma_diag_t magma_diag   = magma_diag_const(diag);

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dtrsm(MagmaLeft, magma_uplo, magma_trans, magma_diag, dim_, dim_, 1.,
        device_data_.get(), ld_, device_data_.get(), ld_,
        magma_singleton.queue_);
}

// get max in absolute value of column j
int ReplicatedMatrix::iamax(const int j, double& val)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    int indx = magma_idamax(dim_, device_data_.get() + j * ld_, 1,
                   magma_singleton.queue_)
               - 1;
    magma_dgetvector(dim_, device_data_.get() + j * ld_ + indx, 1, &val, 1,
        magma_singleton.queue_);

    return indx;
}

void ReplicatedMatrix::setVal(const int i, const int j, const double val)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dsetvector(dim_, &val, 1, device_data_.get() + j * ld_ + i, 1,
        magma_singleton.queue_);
}

void ReplicatedMatrix::setDiagonal(const std::vector<double>& diag_values)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dsetvector(dim_, diag_values.data(), 1, device_data_.get(), ld_ + 1,
        magma_singleton.queue_);
}

double ReplicatedMatrix::norm(char ty)
{
    magma_norm_t magma_ty = magma_norm_const(ty);

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    int lwork = dim_;
    double* dwork;
    magma_dmalloc(&dwork, lwork);
    double norm_val = magmablas_dlange(magma_ty, dim_, dim_, device_data_.get(),
        ld_, dwork, lwork, magma_singleton.queue_);

    magma_singleton.sync();
    magma_free(dwork);

    return norm_val;
}

void ReplicatedMatrix::trset(const char uplo)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    std::vector<double> mat(dim_ * dim_);

    magma_dgetmatrix(dim_, dim_, device_data_.get(), ld_, mat.data(), dim_,
        magma_singleton.queue_);

    if (uplo == 'l' || uplo == 'L')
    {
        for (int j = 0; j < dim_; j++)
            for (int i = 0; i < j; i++)
                mat[i + j * dim_] = 0.;
    }
    else
    {
        for (int j = 0; j < dim_; j++)
            for (int i = j + 1; i < dim_; i++)
                mat[i + j * dim_] = 0.;
    }

    magma_dsetmatrix(dim_, dim_, mat.data(), dim_, device_data_.get(), ld_,
        magma_singleton.queue_);
}

void ReplicatedMatrix::clear()
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magmablas_dlaset(MagmaFull, dim_, dim_, 0.0, 0.0, device_data_.get(), ld_,
        magma_singleton.queue_);
}

void ReplicatedMatrix::print(std::ostream& os, const int ia, const int ja,
    const int ma, const int na) const
{
    const int m = std::min(ma, std::max(dim_ - ia, 0));
    const int n = std::min(na, std::max(dim_ - ja, 0));

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    std::vector<double> mat(dim_ * dim_);

    magma_dgetmatrix(dim_, dim_, device_data_.get(), ld_, mat.data(), dim_,
        magma_singleton.queue_);

    for (int i = ia; i < m; i++)
    {
        for (int j = ja; j < n; j++)
            os << mat[i + j * dim_] << "   ";
        os << std::endl;
    }
}

void ReplicatedMatrix::printMM(std::ostream& os) const {}
#endif
