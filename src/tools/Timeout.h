// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// Adapted form Jeep: Timeout.h,v 1.4 2001/10/10 22:42:44

#ifndef MGMOL_TIMEOUT_H
#define MGMOL_TIMEOUT_H

#include <iomanip>
#include <iostream>
#include <stdlib.h>

#include "MPIdata.h"
#include "mgmol_Signal.h"

#if PCS
#include <csignal>
#include <ctime>
#include <pcserrno.h>
extern "C" int pcssig_register(int signal, time_t mintime, int* status);
void sighandler(int sig);
#endif

#include "MPIdata.h"

class Timeout
{
    int val_;

public:
    int set(const int val)
    {
        val_ = val;
#if PCS
        {
            if (val >= 0)
            {
                if (onpe0)
                {
                    // Register signal on master node only
                    if (val > 0)
                    {
                        int status, stat = 0;
                        stat = pcssig_register(SIGUSR1, val, &status);
                        if (stat != 0)
                        {
                            (*MPIdata::sout)
                                << "TimeoutAlarm: could not register pcs signal"
                                << std::endl;
                            (*MPIdata::sout)
                                << "status = " << status << std::endl;
                        }
                    }
                }
                return 0;
            }
            else
            {
                if (onpe0)
                    (*MPIdata::sout)
                        << " timeout must be non-negative" << std::endl;
                return 1;
            }
        }
#else
        if (onpe0) (*MPIdata::sout) << " timeout not implemented" << std::endl;
        return 1;
#endif
    };

    int value() const { return val_; };

    bool check() { return Signal::received(SIGUSR1); }

    Timeout()
    {
        val_ = 0;
        Signal::enable(SIGUSR1);
    };
};
#endif
