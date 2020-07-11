// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
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

    // blocks CPU until all operations in queue_ are finished
    void sync();

    magma_device_t device_;
    magma_queue_t queue_;

private:
    MagmaSingleton();
};

#endif

#endif
