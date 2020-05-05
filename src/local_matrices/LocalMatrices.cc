// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "LocalMatrices.h"
#include "blas3_c.h"
#include "mputils.h"

template <class T>
LocalMatrices<T>::LocalMatrices(const short subdiv, const int m, const int n)
    : m_(m), n_(n), subdiv_(subdiv)
{
    assert(m >= 0);
    assert(n >= 0);
    assert(subdiv >= 1);

    ptr_matrices_.resize(subdiv_);
    storage_size_ = subdiv_ * m_ * n_;

    if (storage_size_ > 0)
    {
        storage_ = new T[storage_size_];
        memset(storage_, 0, storage_size_ * sizeof(T));
        const int dim2 = m_ * n_;
        for (int i = 0; i < subdiv; i++)
            ptr_matrices_[i] = storage_ + dim2 * i;
    }
    else
    {
        storage_ = nullptr;
    }
}

template <class T>
LocalMatrices<T>::LocalMatrices(const LocalMatrices& mat)
    : m_(mat.m_), n_(mat.n_), subdiv_(mat.subdiv_)
{
    storage_size_ = mat.storage_size_;

    ptr_matrices_.resize(subdiv_);

    if (storage_size_ > 0)
    {
        storage_ = new T[storage_size_];
        memcpy(storage_, mat.storage_, storage_size_ * sizeof(T));
        const int dim2 = m_ * n_;
        for (int i = 0; i < subdiv_; i++)
            ptr_matrices_[i] = storage_ + dim2 * i;
    }
    else
    {
        storage_ = nullptr;
    }
}

template <class T1>
template <class T2>
void LocalMatrices<T1>::copy(const LocalMatrices<T2>& mat)
{
    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        T1* dst_mat             = ptr_matrices_[iloc];
        const T2* const src_mat = mat.getSubMatrix(iloc);
        const int nelements     = n_ * m_;
        for (int j = 0; j < nelements; j++)
        {
            dst_mat[j] = src_mat[j];
        }
    }
}

#ifdef HAVE_BML
template <class T>
void LocalMatrices<T>::copy(const bml_matrix_t* A)
{
    assert(subdiv_ == 1);

    T* dst_mat = ptr_matrices_[0];

    if (bml_get_precision(A) == double_real)
    {
        for (int j = 0; j < n_; j++)
            for (int i = 0; i < m_; i++)
            {
                double* val         = (double*)bml_get(A, i, j);
                dst_mat[i + j * m_] = (T)*val;
            }
    }
    else
    {
        for (int j = 0; j < n_; j++)
            for (int i = 0; i < m_; i++)
            {
                float* val          = (float*)bml_get(A, i, j);
                dst_mat[i + j * m_] = (T)*val;
            }
    }
}
#endif

// perform the symmetric operation
// C := A'*A
// This one is for debug purposes, and can be removed if unused
template <class T>
void LocalMatrices<T>::syrk(
    const int iloc, const int k, const double* const a, const int lda)
{
    assert(iloc < subdiv_);
    assert(iloc < static_cast<int>(ptr_matrices_.size()));
    assert(k <= lda);
    assert(ptr_matrices_[iloc] != nullptr);
    assert(m_ > 0);
    assert(m_ == n_);

    T* const ssiloc = ptr_matrices_[iloc];
    assert(ssiloc != nullptr);

    const double zero = 0.;
    const double one  = 1.;
    const char uplo   = 'l'; // fill lower triangular part
    const char trans  = 't';

    MPsyrk(uplo, trans, m_, k, one, a, lda, zero, ssiloc, m_);
}

// perform the symmetric operation
// C := A'*A
template <class T>
void LocalMatrices<T>::syrk(
    const int iloc, const int k, const float* const a, const int lda)
{
    assert(iloc < subdiv_);
    assert(iloc < (int)ptr_matrices_.size());
    assert(k <= lda);
    assert(ptr_matrices_[iloc] != nullptr);
    assert(m_ > 0);
    assert(m_ == n_);

    T* const ssiloc = ptr_matrices_[iloc];
    assert(ssiloc != nullptr);

    const double zero = 0.;
    const double one  = 1.;
    const char uplo   = 'l'; // fill lower triangular part
    const char trans  = 't';
    MPsyrk(uplo, trans, m_, k, one, a, lda, zero, ssiloc, m_);
    //    dsyrk(&uplo, &trans, &m_, &k, &one, a, &lda,
    //          &zero, ssiloc, &m_);
}

// This one is for debug purposes, and can be removed if unused
template <class T>
void LocalMatrices<T>::gemm(const int iloc, const int ma, const double* const a,
    const int lda, const double* const b, const int ldb)
{
    assert(iloc < subdiv_);
    assert(iloc < (int)ptr_matrices_.size());
    assert(ma <= lda);
    assert(ptr_matrices_[iloc] != nullptr);

    T* const c = ptr_matrices_[iloc];
    assert(c != nullptr);

    MemorySpace::assert_is_host_ptr(a);
    MemorySpace::assert_is_host_ptr(b);
    MemorySpace::assert_is_host_ptr(c);
    LinearAlgebraUtils<MemorySpace::Host>::MPgemm(
        't', 'n', m_, n_, ma, 1., a, lda, b, ldb, 0., c, m_);
}

template <class T>
void LocalMatrices<T>::gemm(const int iloc, const int ma, const float* const a,
    const int lda, const float* const b, const int ldb)
{
    assert(iloc < subdiv_);
    assert(iloc < (int)ptr_matrices_.size());
    assert(ma <= lda);
    assert(ptr_matrices_[iloc] != nullptr);

    T* const c = ptr_matrices_[iloc];
    assert(c != nullptr);

    MemorySpace::assert_is_host_ptr(a);
    MemorySpace::assert_is_host_ptr(b);
    MemorySpace::assert_is_host_ptr(c);
    LinearAlgebraUtils<MemorySpace::Host>::MPgemm(
        't', 'n', m_, n_, ma, 1., a, lda, b, ldb, 0., c, m_);
}

// matrix multiplication
// this = alpha*op(A)*op(B)+beta*this
template <class T>
void LocalMatrices<T>::gemm(const char transa, const char transb,
    const double alpha, const LocalMatrices& matA, const LocalMatrices& matB,
    const double beta)
{
    assert(matA.n() == matB.m());

    const int lda = matA.m();
    const int ldb = matB.m();
    const int nca = matA.n();

    // loop over subdomains
    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        // get matrix data from matA and MatB
        const T* const amat = matA.getSubMatrix(iloc);
        const T* const bmat = matB.getSubMatrix(iloc);
        // get pointer to local storage
        T* const c = ptr_matrices_[iloc];
        assert(c != nullptr);

        // do matrix multiplication
        MemorySpace::assert_is_host_ptr(amat);
        MemorySpace::assert_is_host_ptr(bmat);
        MemorySpace::assert_is_host_ptr(c);
        LinearAlgebraUtils<MemorySpace::Host>::MPgemm(transa, transb, m_, n_,
            nca, alpha, amat, lda, bmat, ldb, beta, c, m_);
    }
}

template <class T>
void LocalMatrices<T>::applyMask(const LocalMatrices& mask)
{
    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        T* local_mat  = ptr_matrices_[iloc];
        T* local_mask = mask.ptr_matrices_[iloc];
        for (int j = 0; j < n_; j++)
        {
            for (int i = 0; i < m_; i++)
            {
                int index = i + j * m_;
                local_mat[index] *= local_mask[index];
            }
        }
    }
}

template <class T>
void LocalMatrices<T>::setMaskThreshold(
    const T min_threshold, const T max_threshold)
{
    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        T* local_mat = ptr_matrices_[iloc];
        for (int j = 0; j < n_; j++)
        {
            for (int i = 0; i < m_; i++)
            {
                int index = i + j * m_;
                if (local_mat[index] > max_threshold
                    || local_mat[index] < min_threshold)
                    local_mat[index] = 0.;
                else
                    local_mat[index] = 1.;
            }
        }
    }
}

template <class T>
void LocalMatrices<T>::printBlock(std::ostream& os, const int blocksize)
{
    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        T* local_mat = ptr_matrices_[iloc];
        os << "iloc=" << iloc << std::endl;
        for (int j = 0; j < blocksize; j++)
        {
            for (int i = 0; i < blocksize; i++)
            {
                int index = i + j * m_;
                os << local_mat[index] << " ";
            }
            os << std::endl;
        }
    }
}

template <class T>
void LocalMatrices<T>::matvec(
    const LocalVector<T>& u, LocalVector<T>& f, const int iloc)
{
    T* mat = ptr_matrices_[iloc];
    Tgemv('n', m_, n_, 1., mat, m_, u.data(), 1, 0., f.data(), 1);
}

template class LocalMatrices<double>;
template class LocalMatrices<float>;

template void LocalMatrices<double>::copy(const LocalMatrices<float>& mat);
template void LocalMatrices<double>::copy(const LocalMatrices<double>& mat);
