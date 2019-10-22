// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_MAGMA_SINGLETON_H
#define MGMOL_MAGMA_SINGLETON_H

#ifdef HAVE_MAGMA
#include <magma_v2.h>

class MagmaSingleton
{
public:
    MagmaSingleton(MagmaSingleton const&) = delete;

    void operator=(MagmaSingleton const&) = delete;

    static MagmaSingleton& get_magma_singleton();

    // This is the real destructor. This should be called before we finalize
    // magma.
    void free();

    magma_device_t device;
    magma_queue_t queue;

private:
    MagmaSingleton();
};

#endif

#endif
