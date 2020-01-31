// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#if defined(HAVE_MAGMA) && defined(HAVE_OPENMP_OFFLOAD)
#include "magma_v2.h"

#include <boost/test/unit_test.hpp>
#include <omp.h>

#include <vector>

BOOST_AUTO_TEST_CASE(magma_openmp)
{
    int on_the_host = -1;
#pragma omp target map(tofrom : on_the_host)
    on_the_host = omp_is_initial_device();
    BOOST_REQUIRE(on_the_host == 0);

    int* val_dev;
    const unsigned int size = 1000;
    magma_malloc(reinterpret_cast<void**>(&val_dev), size * sizeof(int));
    std::vector<int> val_host(size, 1);

    magma_device_t device;
    magma_queue_t queue;
    magma_getdevice(&device);
    magma_queue_create(device, &queue);

    // Copy the values to the device
    magma_setvector(size, sizeof(int), val_host.data(), 1, val_dev, 1, queue);

// Use OpenMP to change the values
#pragma omp target teams distribute parallel for is_device_ptr(val_dev)
    for (unsigned int i = 0; i < size; ++i)
    {
        if (!omp_is_initial_device()) val_dev[i] += i;
    }

    // Copy the values back to the host
    magma_getvector(size, sizeof(int), val_dev, 1, val_host.data(), 1, queue);

    // Check the result
    for (unsigned int i = 0; i < size; ++i)
        BOOST_TEST(val_host[i] == i + 1);

    magma_queue_destroy(queue);
}
#endif
