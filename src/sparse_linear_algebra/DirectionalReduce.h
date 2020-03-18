// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include <iostream>
#include <mpi.h>

class DirectionalReduce
{
    const MPI_Comm cart_comm_;

    MPI_Request request_[4];

    int sbuf_[2]; // work buffer for data transfer (send)
    int rbuf_[2]; // work buffer for data transfer (recv)

    int lstep_[3];
    int rstep_[3];

    int nprocs_[3];

    void setupPersistentRequests(const short dir);

    void deletePersistentRequests();

    void computeNumSteps(const int max_steps[3]);
    void computeNumSteps(const double spread_radius, const double domain[3]);

public:
    DirectionalReduce(MPI_Comm cart_comm, const int max_steps[3]);

    DirectionalReduce(
        MPI_Comm cart_comm, const double spread_radius, const double domain[3]);

    void computeDirMax(const short dir, int data[2]);

    int lstep(const short dir) { return lstep_[dir]; }
    int rstep(const short dir) { return rstep_[dir]; }

    void printStats(std::ostream& os);
};
