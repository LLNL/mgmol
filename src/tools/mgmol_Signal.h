// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// Adapted from Jeep: Signal.h,v 1.5 2002/06/28 20:50:33

// The Signal class is a utility to catch UNIX signals.
// A set of flags is maintained to remeber which signals were caught.
// Member functions are provided to enable or disable signals,
// access the flag set, reset flags, or interrogate flags.
// The Signal class can be used in an application by declaring
//
// #include "mgmol_Signal.h"
// set<int> Signal::recv_;
//
// A signal can be registered using, e.g.
//
// Signal::enable(SIGUSR1);
//
// To test if signal SIGUSR1 has been received, use
// if ( Signal::received(SIGUSR1) ) { ... }

#ifndef MGMOL_SIGNAL_H
#define MGMOL_SIGNAL_H

#include <csignal>
#include <mpi.h>
#include <set>

class Signal
{
    static std::set<int> recv_;

public:
    static void enable(int sig) { signal(sig, (void (*)(int))(&setsig)); }
    static void disable(int sig)
    {
        signal(sig, SIG_IGN);
        reset(sig);
    }
    static void setsig(int sig) { recv_.insert(sig); }
    static void reset(int sig) { recv_.erase(sig); }

    static bool received(int sig)
    {
        int flag = recv_.find(sig) != recv_.end();

        MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

        return flag != 0;
    }
};

#endif
