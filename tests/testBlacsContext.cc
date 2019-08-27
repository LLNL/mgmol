// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "../DistMatrix/BlacsContext.h"

#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_AUTO_TEST_CASE(blacs)
{
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // build a BlacsContext 1 row x nprocs col
    dist_matrix::BlacsContext bc(MPI_COMM_WORLD, 'r', nprocs);

    BOOST_TEST(bc.nprow() == 1);
    BOOST_TEST(bc.npcol() == nprocs);
    BOOST_TEST(bc.nprocs() == nprocs);
    BOOST_TEST(bc.myproc() < nprocs);
    BOOST_TEST(bc.myproc() >= 0);

    int sendbuff = bc.myproc();
    int recvbuf  = 0;
    int sum      = nprocs * (nprocs + 1) / 2 - nprocs;
    MPI_Allreduce(&sendbuff, &recvbuf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    BOOST_TEST(recvbuf == sum);
}
