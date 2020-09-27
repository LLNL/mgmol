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

#include "magma_v2.h"

using MemoryDev = MemorySpace::Memory<double, MemorySpace::Device>;

ReplicatedVector::ReplicatedVector(const std::string name, const int n)
    : dim_(n),
      device_data_(MemoryDev::allocate(dim_), MemoryDev::free)
{
}

ReplicatedVector::ReplicatedVector(const ReplicatedVector& v)
    : dim_(v.dim_),
      device_data_(MemoryDev::allocate(dim_), MemoryDev::free)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dcopy(dim_, v.device_data_.get(), 1, device_data_.get(), 1,
        magma_singleton.queue_);
}

ReplicatedVector::ReplicatedVector(const std::vector<double>& v)
    : dim_(v.size()),
      device_data_(MemoryDev::allocate(dim_), MemoryDev::free)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dsetvector(dim_, v.data(), 1, device_data_.get(), 1, magma_singleton.queue_);
}

ReplicatedVector& ReplicatedVector::operator=(const ReplicatedVector& src)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_dcopy(dim_, src.device_data_.get(), 1, device_data_.get(), 1,
        magma_singleton.queue_);

    return *this;
}

void ReplicatedVector::axpy(const double alpha, const ReplicatedVector& x)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magma_daxpy(dim_, alpha, x.device_data_.get(), 1, device_data_.get(), 1,
        magma_singleton.queue_);
}

void ReplicatedVector::clear()
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    std::vector<double> zero(dim_,0.);

    magma_dsetvector(dim_, zero.data(), 1, device_data_.get(), 1, magma_singleton.queue_);
}

double ReplicatedVector::nrm2()
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    return magma_dnrm2(dim_, device_data_.get(), 1, magma_singleton.queue_);
}

double ReplicatedVector::dot(const ReplicatedVector& v)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    return magma_ddot(dim_, device_data_.get(), 1, 
                      v.device_data_.get(), 1, magma_singleton.queue_);
}

void ReplicatedVector::gemv(const char trans,
    const double alpha, const ReplicatedMatrix& a, const ReplicatedVector& b,
    const double beta)
{
    magma_trans_t magma_trans = magma_trans_const(trans);

    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    magmablas_dgemv(magma_trans, dim_, dim_, alpha,
        a.device_data_.get(), a.ld_, b.device_data_.get(), 1, beta,
        device_data_.get(), 1, magma_singleton.queue_);
}

