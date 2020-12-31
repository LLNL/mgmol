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
#include "magma_singleton.h"
#include "mputils.h"

template <typename DataType, typename MemorySpaceType>
LocalMatrices<DataType, MemorySpaceType>::LocalMatrices(
    const short nmat, const int m, const int n)
    : m_(m),
      n_(n),
      nmat_(nmat),
      storage_(MemorySpace::Memory<DataType, MemorySpaceType>::allocate(
                   nmat_ * m_ * n_),
          MemorySpace::Memory<DataType, MemorySpaceType>::free)
{
    assert(m >= 0);
    assert(n >= 0);
    assert(nmat >= 1);

    allocate();

    if (storage_size_ > 0)
    {
        set2zero();
    }
}

template <>
void LocalMatrices<double, MemorySpace::Host>::set2zero()
{
    memset(storage_.get(), 0, storage_size_ * sizeof(double));
}

template <>
void LocalMatrices<float, MemorySpace::Host>::set2zero()
{
    memset(storage_.get(), 0, storage_size_ * sizeof(float));
}

#ifdef HAVE_MAGMA
template <>
void LocalMatrices<double, MemorySpace::Device>::set2zero()
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magmablas_dlaset(MagmaFull, m_, n_ * nmat_, 0.0, 0.0, storage_.get(), m_,
        magma_singleton.queue_);
}

template <>
void LocalMatrices<float, MemorySpace::Device>::set2zero()
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magmablas_slaset(MagmaFull, m_, n_ * nmat_, 0.0, 0.0, storage_.get(), m_,
        magma_singleton.queue_);
}
#endif

template <typename DataType, typename MemorySpaceType>
LocalMatrices<DataType, MemorySpaceType>::LocalMatrices(
    const LocalMatrices& mat)
    : m_(mat.m_),
      n_(mat.n_),
      nmat_(mat.nmat_),
      storage_(MemorySpace::Memory<DataType, MemorySpaceType>::allocate(
                   nmat_ * m_ * n_),
          MemorySpace::Memory<DataType, MemorySpaceType>::free)
{
    allocate();

    if (storage_size_ > 0)
    {
        memcpy(storage_.get(), mat.storage_.get(),
            storage_size_ * sizeof(DataType));
    }
}

#ifdef HAVE_MAGMA
template <>
LocalMatrices<double, MemorySpace::Device>::LocalMatrices(
    const LocalMatrices& mat)
    : m_(mat.m_),
      n_(mat.n_),
      nmat_(mat.nmat_),
      storage_(MemorySpace::Memory<double, MemorySpace::Device>::allocate(
                   nmat_ * m_ * n_),
          MemorySpace::Memory<double, MemorySpace::Device>::free)
{
    allocate();

    if (storage_size_ > 0)
    {
        auto& magma_singleton = MagmaSingleton::get_magma_singleton();
        magma_dcopymatrix(m_, n_, mat.storage_.get(), m_, storage_.get(), m_,
            magma_singleton.queue_);
    }
}

template <>
LocalMatrices<float, MemorySpace::Device>::LocalMatrices(
    const LocalMatrices& mat)
    : m_(mat.m_),
      n_(mat.n_),
      nmat_(mat.nmat_),
      storage_(MemorySpace::Memory<float, MemorySpace::Device>::allocate(
                   nmat_ * m_ * n_),
          MemorySpace::Memory<float, MemorySpace::Device>::free)
{
    allocate();

    if (storage_size_ > 0)
    {
        auto& magma_singleton = MagmaSingleton::get_magma_singleton();
        magma_scopymatrix(m_, n_, mat.storage_.get(), m_, storage_.get(), m_,
            magma_singleton.queue_);
    }
}
#endif

template <typename DataType, typename MemorySpaceType>
void LocalMatrices<DataType, MemorySpaceType>::allocate()
{
    ptr_matrices_.resize(nmat_);
    storage_size_ = nmat_ * m_ * n_;

    if (storage_size_ > 0)
    {
        const int dim2 = m_ * n_;
        for (int i = 0; i < nmat_; i++)
            ptr_matrices_[i] = storage_.get() + dim2 * i;
    }
}

template <typename DataType, typename MemorySpaceType>
template <typename DataType2>
void LocalMatrices<DataType, MemorySpaceType>::copy(
    const LocalMatrices<DataType2, MemorySpaceType>& mat)
{
    for (short iloc = 0; iloc < nmat_; iloc++)
    {
        DataType* dst_mat              = ptr_matrices_[iloc];
        const DataType2* const src_mat = mat.getSubMatrix(iloc);
        const int nelements            = n_ * m_;
        for (int j = 0; j < nelements; j++)
        {
            dst_mat[j] = src_mat[j];
        }
    }
}

#ifdef HAVE_BML
template <typename DataType, typename MemorySpaceType>
void LocalMatrices<DataType, MemorySpaceType>::copy(const bml_matrix_t* A)
{
    assert(nmat_ == 1);

    DataType* dst_mat = ptr_matrices_[0];

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

template <>
void LocalMatrices<double, MemorySpace::Host>::scal(const double alpha)
{
    LinearAlgebraUtils<MemorySpace::Host>::MPscal(
        storage_size_, alpha, storage_.get());
}

#ifdef HAVE_MAGMA
template <>
void LocalMatrices<double, MemorySpace::Device>::scal(const double alpha)
{
    LinearAlgebraUtils<MemorySpace::Device>::MPscal(
        storage_size_, alpha, storage_.get());
}
#endif

// perform the symmetric operation
// C := A'*A
// This one is for debug purposes, and can be removed if unused
template <>
void LocalMatrices<double, MemorySpace::Host>::syrk(
    const int iloc, const int k, const double* const a, const int lda)
{
    assert(iloc < nmat_);
    assert(iloc < static_cast<int>(ptr_matrices_.size()));
    assert(k <= lda);
    assert(ptr_matrices_[iloc] != nullptr);
    assert(m_ > 0);
    assert(m_ == n_);

    double* const ssiloc = ptr_matrices_[iloc];
    assert(ssiloc != nullptr);

    const double zero = 0.;
    const double one  = 1.;
    const char uplo   = 'l'; // fill lower triangular part
    const char trans  = 't';

    // ssiloc is always on the host but a maybe on the device
#ifdef HAVE_MAGMA
    const int size_a = lda * m_;
    double* a_host_view
        = MemorySpace::Memory<double, MemorySpace::Device>::allocate_host_view(
            size_a);
    MemorySpace::Memory<double, MemorySpace::Device>::copy_view_to_host(
        const_cast<double*>(a), size_a, a_host_view);
#else
    const double* const a_host_view = a;
#endif

    LinearAlgebraUtils<MemorySpace::Host>::MPsyrk(
        uplo, trans, m_, k, one, a_host_view, lda, zero, ssiloc, m_);

#ifdef HAVE_MAGMA
    MemorySpace::Memory<double, MemorySpace::Device>::free_host_view(
        a_host_view);
#endif
}

#ifdef HAVE_MAGMA
template <>
void LocalMatrices<double, MemorySpace::Device>::syrk(
    const int iloc, const int k, const double* const a, const int lda)
{
    assert(iloc < nmat_);
    assert(iloc < static_cast<int>(ptr_matrices_.size()));
    assert(k <= lda);
    assert(ptr_matrices_[iloc] != nullptr);
    assert(m_ > 0);
    assert(m_ == n_);
    MemorySpace::assert_is_host_ptr(a);

    double* const ssiloc = ptr_matrices_[iloc];
    assert(ssiloc != nullptr);

    const double zero = 0.;
    const double one  = 1.;
    const char uplo   = 'l'; // fill lower triangular part
    const char trans  = 't';

    LinearAlgebraUtils<MemorySpace::Device>::MPsyrk(
        uplo, trans, m_, k, one, a, lda, zero, ssiloc, m_);
}
#endif

// perform the symmetric operation
// C := A'*A
template <typename DataType, typename MemorySpaceType>
void LocalMatrices<DataType, MemorySpaceType>::syrk(
    const int iloc, const int k, const float* const a, const int lda)
{
    assert(iloc < nmat_);
    assert(iloc < (int)ptr_matrices_.size());
    assert(k <= lda);
    assert(ptr_matrices_[iloc] != nullptr);
    assert(m_ > 0);
    assert(m_ == n_);

    DataType* const ssiloc = ptr_matrices_[iloc];
    assert(ssiloc != nullptr);

    const double zero = 0.;
    const double one  = 1.;
    const char uplo   = 'l'; // fill lower triangular part
    const char trans  = 't';
    LinearAlgebraUtils<MemorySpaceType>::MPsyrk(
        uplo, trans, m_, k, one, a, lda, zero, ssiloc, m_);
}

// This one is for debug purposes, and can be removed if unused
template <typename DataType, typename MemorySpaceType>
void LocalMatrices<DataType, MemorySpaceType>::gemm(const int iloc,
    const int ma, const double* const a, const int lda, const double* const b,
    const int ldb)
{
    assert(iloc < nmat_);
    assert(iloc < (int)ptr_matrices_.size());
    assert(ma <= lda);
    assert(ptr_matrices_[iloc] != nullptr);

    DataType* const c = ptr_matrices_[iloc];
    assert(c != nullptr);

    MemorySpace::assert_is_host_ptr(a);
    MemorySpace::assert_is_host_ptr(b);
    MemorySpace::assert_is_host_ptr(c);
    LinearAlgebraUtils<MemorySpaceType>::MPgemm(
        't', 'n', m_, n_, ma, 1., a, lda, b, ldb, 0., c, m_);
}

template <typename DataType, typename MemorySpaceType>
void LocalMatrices<DataType, MemorySpaceType>::gemm(const int iloc,
    const int ma, const float* const a, const int lda, const float* const b,
    const int ldb)
{
    assert(iloc < nmat_);
    assert(iloc < (int)ptr_matrices_.size());
    assert(ma <= lda);
    assert(ptr_matrices_[iloc] != nullptr);

    DataType* const c = ptr_matrices_[iloc];
    assert(c != nullptr);

    MemorySpace::assert_is_host_ptr(a);
    MemorySpace::assert_is_host_ptr(b);
    MemorySpace::assert_is_host_ptr(c);
    LinearAlgebraUtils<MemorySpaceType>::MPgemm(
        't', 'n', m_, n_, ma, 1., a, lda, b, ldb, 0., c, m_);
}

// matrix multiplication
// this = alpha*op(A)*op(B)+beta*this
template <typename DataType, typename MemorySpaceType>
void LocalMatrices<DataType, MemorySpaceType>::gemm(const char transa,
    const char transb, const double alpha, const LocalMatrices& matA,
    const LocalMatrices& matB, const double beta)
{
    assert(matA.n() == matB.m());

    const int lda = matA.m();
    const int ldb = matB.m();
    const int nca = matA.n();

    // loop over subdomains
    for (short iloc = 0; iloc < nmat_; iloc++)
    {
        // get matrix data from matA and MatB
        const DataType* const amat = matA.getSubMatrix(iloc);
        const DataType* const bmat = matB.getSubMatrix(iloc);
        // get pointer to local storage
        DataType* const c = ptr_matrices_[iloc];
        assert(c != nullptr);

        // do matrix multiplication
        MemorySpace::assert_is_host_ptr(amat);
        MemorySpace::assert_is_host_ptr(bmat);
        MemorySpace::assert_is_host_ptr(c);
        LinearAlgebraUtils<MemorySpaceType>::MPgemm(transa, transb, m_, n_, nca,
            alpha, amat, lda, bmat, ldb, beta, c, m_);
    }
}

template <typename DataType, typename MemorySpaceType>
void LocalMatrices<DataType, MemorySpaceType>::applyMask(
    const LocalMatrices& mask)
{
    for (short iloc = 0; iloc < nmat_; iloc++)
    {
        DataType* local_mat  = ptr_matrices_[iloc];
        DataType* local_mask = mask.ptr_matrices_[iloc];
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

template <typename DataType, typename MemorySpaceType>
void LocalMatrices<DataType, MemorySpaceType>::setMaskThreshold(
    const DataType min_threshold, const DataType max_threshold)
{
    for (short iloc = 0; iloc < nmat_; iloc++)
    {
        DataType* local_mat = ptr_matrices_[iloc];
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

template <>
void LocalMatrices<double, MemorySpace::Host>::setValues(
    double* values, const int ld, const int iloc)
{
    double* local_mat = ptr_matrices_[iloc];
    for (int j = 0; j < n_; j++)
    {
        memcpy(local_mat + j * m_, values + j * ld, n_ * sizeof(double));
    }
}

#ifdef HAVE_MAGMA
template <>
void LocalMatrices<double, MemorySpace::Device>::setValues(
    double* values, const int ld, const int iloc)
{
    MemorySpace::assert_is_dev_ptr(values);

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();
    magma_dcopymatrix(
        m_, n_, values, ld, ptr_matrices_[iloc], m_, magma_singleton.queue_);
}

template <>
void LocalMatrices<double, MemorySpace::Device>::assign(
    const LocalMatrices<double, MemorySpace::Host>& src)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    for (short iloc = 0; iloc < nmat_; iloc++)
        magma_dsetmatrix(src.m(), src.n(), src.getSubMatrix(iloc), src.n(),
            ptr_matrices_[iloc], m_, magma_singleton.queue_);
}

template <>
void LocalMatrices<double, MemorySpace::Host>::assign(
    const LocalMatrices<double, MemorySpace::Device>& src)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    for (short iloc = 0; iloc < nmat_; iloc++)
        magma_dgetmatrix(src.m(), src.n(), src.getSubMatrix(iloc), src.n(),
            ptr_matrices_[iloc], m_, magma_singleton.queue_);
}
#endif

template <typename DataType, typename MemorySpaceType>
void LocalMatrices<DataType, MemorySpaceType>::printBlock(
    std::ostream& os, const int blocksize)
{
    for (short iloc = 0; iloc < nmat_; iloc++)
    {
        DataType* local_mat = ptr_matrices_[iloc];
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

template <typename DataType, typename MemorySpaceType>
void LocalMatrices<DataType, MemorySpaceType>::matvec(
    const LocalVector<DataType, MemorySpaceType>& u,
    LocalVector<DataType, MemorySpaceType>& f, const int iloc)
{
    DataType* mat = ptr_matrices_[iloc];
    Tgemv('n', m_, n_, 1., mat, m_, u.data(), 1, 0., f.data(), 1);
}

template class LocalMatrices<double, MemorySpace::Host>;
template class LocalMatrices<float, MemorySpace::Host>;

template void LocalMatrices<double, MemorySpace::Host>::copy(
    const LocalMatrices<float, MemorySpace::Host>& mat);
template void LocalMatrices<double, MemorySpace::Host>::copy(
    const LocalMatrices<double, MemorySpace::Host>& mat);

#ifdef HAVE_MAGMA
template class LocalMatrices<double, MemorySpace::Device>;
#endif
