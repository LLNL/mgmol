// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "../DistMatrix/BlacsContext.h"

#include "catch.hpp"
#include <iostream>

TEST_CASE("Simple check for blacs", "[blacs]")
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // build a BlacsContext 1 row x nprocs col
    dist_matrix::BlacsContext bc(MPI_COMM_WORLD, 'r', nprocs);

    CHECK(bc.nprow() == 1);
    CHECK(bc.npcol() == nprocs);
    CHECK(bc.nprocs() == nprocs);
    CHECK(bc.myproc() < nprocs);
    CHECK(bc.myproc() >= 0);

    int sendbuff = bc.myproc();
    int recvbuf  = 0;
    int sum      = nprocs * (nprocs + 1) / 2 - nprocs;
    MPI_Allreduce(&sendbuff, &recvbuf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CHECK(recvbuf == sum);
}
