// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
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

#include <cstring>
#include <memory>
#include <vector>

#define MGMOL_STRINGIFY_(...) #__VA_ARGS__
#define MGMOL_STRINGIFY(...) MGMOL_STRINGIFY_(__VA_ARGS__)

#if defined(HAVE_MAGMA) && defined(HAVE_OPENMP_OFFLOAD)
#define MGMOL_PARALLEL_FOR(...)                                                \
  _Pragma(MGMOL_STRINGIFY(omp target teams distribute parallel for is_device_ptr(__VA_ARGS__)))
#define MGMOL_PARALLEL_FOR_COLLAPSE(n, ...)                                    \
  _Pragma(MGMOL_STRINGIFY(omp target teams distribute parallel for collapse(n) is_device_ptr(__VA_ARGS__)))
#else
#define MGMOL_PARALLEL_FOR(...) _Pragma(MGMOL_STRINGIFY(omp parallel for))
#define MGMOL_PARALLEL_FOR_COLLAPSE(n, ...) _Pragma(MGMOL_STRINGIFY(omp parallel for collapse (n)))
#endif

namespace MemorySpace
{
struct Host
{
};

struct Device
{
};

//---------------------------------------------------------------------------//
// Copy from the host/device to the device/host
//---------------------------------------------------------------------------//

#ifdef HAVE_MAGMA
template <typename T>
void copy_to_dev(T const* const vec, unsigned int size, T* vec_dev)
{
    int const increment   = 1;
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();
    magma_setvector(size, sizeof(T), vec, increment, vec_dev, increment,
        magma_singleton.queue_);
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
        magma_singleton.queue_);
}
#else
template <typename T>
void copy_to_host(T const* const /*vec_dev*/, unsigned int /*size*/, T* /*vec*/)
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

//---------------------------------------------------------------------------//
// Allocate and deallocate memory
//---------------------------------------------------------------------------//

template <typename MemorySpaceType>
struct Memory
{
    template <typename T>
    static T* allocate(unsigned int size);

    template <typename T>
    static void free(T* ptr);

    template <typename T>
    static void copy(T const* in, unsigned int size, T* out);

    template <typename T>
    static void set(T* ptr, unsigned int size, int val);
};

template <>
struct Memory<MemorySpace::Host>
{
    template <typename T>
    static T* allocate(unsigned int size)
    {
        return std::malloc(size * sizeof(T));
    }

    template <typename T>
    static void free(T* ptr)
    {
        std::free(ptr);
        ptr = nullptr;
    }

    template <typename T>
    static void copy(T const* in, unsigned int size, T* out)
    {
        std::memcpy(out, in, size * sizeof(T));
    }

    template <typename T>
    static void set(T* ptr, unsigned int size, int val)
    {
        std::memset(ptr, val, size * sizeof(T));
    }
};

#ifdef HAVE_MAGMA
template <>
struct Memory<MemorySpace::Device>
{
    template <typename T>
    static T* allocate(unsigned int size)
    {
        T* ptr_dev;
        magma_malloc(reinterpret_cast<void**>(&ptr_dev), size * sizeof(T));

        return ptr_dev;
    }

    template <typename T>
    static void free(T* ptr_dev)
    {
        magma_free(ptr_dev);
        ptr_dev = nullptr;
    }

    template <typename T>
    static void copy(T const* in, unsigned int size, T* out)
    {
        int const increment   = 1;
        auto& magma_singleton = MagmaSingleton::get_magma_singleton();
        magma_copyvector(size, sizeof(T), in, increment, out, increment,
            magma_singleton.queue_);
    }

    template <typename T>
    static void set(T* ptr, unsigned int size, int val)
    {
        // Cannot directly use cudaMemset and MAGMA does not have an equivalent
        // function. If we have OpenMP with offloading, we directly set the
        // value. Otherwise, we copy a vector from the host.
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp target teams distribute parallel for is_device_ptr(ptr)
        for (unsigned int i = 0; i < size; ++i)
            ptr[i] = val;
#else
        auto ptr_host = Memory<Host>::allocate<T>(size);
        set(ptr_host, size, val);
        copy_to_dev(ptr_host, size, ptr);
#endif
    }
};
#endif
}
#endif
