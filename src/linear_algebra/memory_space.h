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

#ifdef HAVE_MAGMA
#include <magma_v2.h>
#endif

#include <magma_singleton.h>

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

#ifdef HAVE_MAGMA
template <typename T>
void copy_to_dev(T const* const vec, unsigned int size, T* vec_dev)
{
    int const increment   = 1;
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();
    magma_setvector(size, sizeof(T), vec, increment, vec_dev, increment,
        magma_singleton.queue);
}
#else

template <typename T>
void copy_to_dev(T* const, unsigned int, T*)
{
    // TODO
}
#endif

template <typename T>
void copy_to_dev(std::vector<T> const& vec, T* vec_dev)
{
    copy_to_dev(vec.data(), vec.size(), vec_dev);
}

template <typename T>
void copy_to_dev(
    std::vector<T> const& vec, std::unique_ptr<T[], void (*)(T*)>& vec_dev)
{
    copy_to_dev(vec.data(), vec.size(), vec_dev.get());
}

template <typename T>
void copy_to_dev(std::vector<T> const& vec, std::shared_ptr<T[]>& vec_dev)
{
    copy_to_dev(vec, vec.size, vec_dev.get());
}

template <typename T>
void copy_to_dev(T const* const vec, unsigned int size,
    std::unique_ptr<T[], void (*)(T*)>& vec_dev)
{
    copy_to_dev(vec, size, vec_dev.get());
}

#ifdef HAVE_MAGMA
template <typename T>
void copy_to_host(T const* const vec_dev, unsigned int size, T* vec)
{
    int const increment   = 1;
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();
    magma_getvector(size, sizeof(T), vec_dev, increment, vec, increment,
        magma_singleton.queue);
}
#else
template <typename T>
void copy_to_host(T const* const vec_dev, unsigned int size, T* vec)
{
    // TODO
}
#endif

template <typename T>
void copy_to_host(T const* const vec_dev, std::vector<T>& vec)
{
    copy_to_host(vec_dev, vec.size(), vec.data());
}

template <typename T>
void copy_to_host(
    std::unique_ptr<T[], void (*)(T*)> const& vec_dev, std::vector<T>& vec)
{
    copy_to_host(vec_dev.get(), vec.size(), vec.data());
}

template <typename T>
void copy_to_host(std::shared_ptr<T[]> const& vec_dev, std::vector<T>& vec)
{
    copy_to_host(vec_dev.get(), vec.size(), vec.data());
}

template <typename T>
void copy_to_host(std::unique_ptr<T[], void (*)(T*)> const& vec_dev,
    unsigned int size, T* vec)
{
    copy_to_host(vec_dev.get(), size, vec);
}

#ifdef HAVE_MAGMA
template <typename T>
T* allocate_data_dev(unsigned int size)
{
    T* ptr_dev;
    magma_malloc(reinterpret_cast<void**>(&ptr_dev), size * sizeof(T));

    return ptr_dev;
}

template <typename T>
void delete_data_dev(T* ptr_dev)
{
    magma_free(ptr_dev);
}
#endif
}
#endif
