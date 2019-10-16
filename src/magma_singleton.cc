// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <magma_singleton.h>

#ifdef HAVE_MAGMA

MagmaSingleton::MagmaSingleton()
{
    magma_getdevice(&device);
    magma_queue_create(device, &queue);
}

MagmaSingleton& MagmaSingleton::get_magma_singleton()
{
    static MagmaSingleton singleton;

    return singleton;
}

void MagmaSingleton::free() { magma_queue_destroy(queue); }

#endif
