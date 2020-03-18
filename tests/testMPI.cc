// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "catch.hpp"
#include <iostream>
#include <mpi.h>

TEST_CASE("MPI sanity check", "[MPI]")
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    CHECK(nprocs == 4);

    int myproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

    int sendbuff = myproc;
    int recvbuf  = 0;
    int sum      = nprocs * (nprocs + 1) / 2 - nprocs;
    int mpi_ret  = MPI_Allreduce(
        &sendbuff, &recvbuf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK(mpi_ret != MPI_ERR_COMM);
    CHECK(recvbuf == sum);
}
