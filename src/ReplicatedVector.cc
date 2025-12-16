// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "ReplicatedVector.h"

#include "memory_space.h"

#ifdef HAVE_MAGMA
#include "magma_v2.h"
#else
#include "blas2_c.h"
#endif

#ifdef HAVE_MAGMA
using MemoryDev = MemorySpace::Memory<double, MemorySpace::Device>;
#else
using MemoryDev = MemorySpace::Memory<double, MemorySpace::Host>;
#endif

ReplicatedVector::ReplicatedVector(const int n)
    : dim_(n), data_(MemoryDev::allocate(dim_), MemoryDev::free)
{
}

ReplicatedVector::ReplicatedVector(const ReplicatedVector& v)
    : dim_(v.dim_), data_(MemoryDev::allocate(dim_), MemoryDev::free)
{
#ifdef HAVE_MAGMA
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dcopy(dim_, v.data_.get(), 1, data_.get(), 1, magma_singleton.queue_);
#else
    memcpy(data_.get(), v.data_.get(), dim_ * sizeof(double));
#endif
}

ReplicatedVector::ReplicatedVector(const std::vector<double>& v)
    : dim_(v.size()), data_(MemoryDev::allocate(dim_), MemoryDev::free)
{
#ifdef HAVE_MAGMA
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dsetvector(dim_, v.data(), 1, data_.get(), 1, magma_singleton.queue_);
#else
    memcpy(data_.get(), v.data(), dim_ * sizeof(double));
#endif
}

ReplicatedVector& ReplicatedVector::operator=(const ReplicatedVector& src)
{
#ifdef HAVE_MAGMA
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dcopy(
        dim_, src.data_.get(), 1, data_.get(), 1, magma_singleton.queue_);
#else
    memcpy(data_.get(), src.data_.get(), dim_ * sizeof(double));
#endif

    return *this;
}

void ReplicatedVector::axpy(const double alpha, const ReplicatedVector& x)
{
#ifdef HAVE_MAGMA
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_daxpy(
        dim_, alpha, x.data_.get(), 1, data_.get(), 1, magma_singleton.queue_);
#else
    int ione = 1;
    double a = alpha;
    DAXPY(&dim_, &a, x.data_.get(), &ione, data_.get(), &ione);
#endif
}

void ReplicatedVector::clear()
{
#ifdef HAVE_MAGMA
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    std::vector<double> zero(dim_, 0.);

    magma_dsetvector(
        dim_, zero.data(), 1, data_.get(), 1, magma_singleton.queue_);
#else
    memset(data_.get(), 0., dim_ * sizeof(double));
#endif
}

double ReplicatedVector::nrm2()
{
#ifdef HAVE_MAGMA
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    return magma_dnrm2(dim_, data_.get(), 1, magma_singleton.queue_);
#else
    int ione = 1;
    return DNRM2(&dim_, data_.get(), &ione);
#endif
}

double ReplicatedVector::dot(const ReplicatedVector& v)
{
#ifdef HAVE_MAGMA
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    return magma_ddot(
        dim_, data_.get(), 1, v.data_.get(), 1, magma_singleton.queue_);
#else
    int ione = 1;
    return DDOT(&dim_, data_.get(), &ione, v.data_.get(), &ione);
#endif
}

void ReplicatedVector::gemv(const char trans, const double alpha,
    const ReplicatedMatrix& a, const ReplicatedVector& b, const double beta)
{
#ifdef HAVE_MAGMA
    magma_trans_t magma_trans = magma_trans_const(trans);

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magmablas_dgemv(magma_trans, dim_, dim_, alpha, a.data_.get(), a.ld_,
        b.data_.get(), 1, beta, data_.get(), 1, magma_singleton.queue_);
#else
    int ione = 1;
    int lda  = a.ld_;
    DGEMV(&trans, &dim_, &dim_, &alpha, a.data_.get(), &lda, b.data_.get(),
        &ione, &beta, data_.get(), &ione);
#endif
}
