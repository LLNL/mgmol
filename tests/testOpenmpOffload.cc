// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "catch.hpp"
#include "memory_space.h"
#include <omp.h>

#include "magma_v2.h"

template <typename ScalarType>
using MemoryDev = MemorySpace::Memory<ScalarType, MemorySpace::Device>;

TEST_CASE("Check OpenMP offload wrapper", "[openmp_wrapper]")
{
    magma_device_t device;
    magma_queue_t queue;
    magma_getdevice(&device);
    magma_queue_create(device, &queue);

    const unsigned int size = 1024;

    std::vector<int> val_host(size, 1);

#ifdef HAVE_OPENMP_OFFLOAD
    std::unique_ptr<int, void (*)(int*)> val_dev(
        MemoryDev<int>::allocate(size), MemoryDev<int>::free);

    MemorySpace::copy_to_dev(val_host, val_dev);

    int* val_alias = val_dev.get();
#else
    int* val_alias = val_host.data();
#endif

    // Use OpenMP to change the values
    MGMOL_PARALLEL_FOR(val_alias)
    for (unsigned int i = 0; i < size; ++i)
    {
        val_alias[i] += i;
    }

#ifdef HAVE_OPENMP_OFFLOAD
    // Copy the values back to the host
    MemorySpace::copy_to_host(val_alias, val_host);
#endif

    // Check the result
    for (unsigned int i = 0; i < size; ++i)
        CHECK(val_host[i] == i + 1);

    magma_queue_destroy(queue);
}
