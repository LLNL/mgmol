// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_MEMORY_SPACE_H
#define MGMOL_MEMORY_SPACE_H

#ifdef USE_MAGMA
#include <magma_v2.h>
#endif

#include <memory>
#include <vector>

namespace MemorySpace
{
struct Host
{
};

struct Device
{
};

#ifdef USE_MAGMA
template <typename T>
void copy_to_dev(std::vector<T> const& vec, T* vec_dev)
{
    // TODO this should be moved to a singleton but we should call the
    // destructor before calling magma finalize
    magma_device_t device;
    magma_queue_t queue;
    magma_getdevice(&device);
    magma_queue_create(device, &queue);

    magma_setvector(vec.size(), sizeof(T), vec.data(), 1, vec_dev, 1, queue);

    // TODO this should be in the singleton
    magma_queue_destroy(queue);
}
#else
template <typename T>
void copy_to_dev(std::vector<T> const&, T*)
{
}
#endif

template <typename T>
void copy_to_dev(
    std::vector<T> const& vec, std::unique_ptr<T[], void (*)(T*)>& vec_dev)
{
    copy_to_dev(vec, vec_dev.get());
}

template <typename T>
void copy_to_dev(std::vector<T> const& vec, std::shared_ptr<T[]>& vec_dev)
{
    copy_to_dev(vec, vec_dev.get());
}

#ifdef USE_MAGMA
template <typename T>
void copy_to_host(T const* const vec_dev, std::vector<T>& vec)
{
    // TODO this should be moved to a singleton but we should call the
    // destructor before calling magma finalize
    magma_device_t device;
    magma_queue_t queue;
    magma_getdevice(&device);
    magma_queue_create(device, &queue);

    magma_getvector(vec.size(), sizeof(T), vec_dev, 1, vec.data(), 1, queue);

    // TODO this should be in the singleton
    magma_queue_destroy(queue);
}
#else
template <typename T>
void copy_to_host(T const* const, std::vector<T>&)
{
}
#endif

template <typename T>
void copy_to_host(
    std::unique_ptr<T[], void (*)(T*)> const& vec_dev, std::vector<T>& vec)
{
    copy_to_host(vec_dev.get(), vec);
}

template <typename T>
void copy_to_host(std::shared_ptr<T[]> const& vec_dev, std::vector<T>& vec)
{
    copy_to_host(vec_dev.get(), vec);
}
}
#endif
