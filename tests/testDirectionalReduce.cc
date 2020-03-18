// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "DirectionalReduce.h"
#include "PEenv.h"

#include "catch.hpp"
#include <iostream>

TEST_CASE("Check DirectionalReduce", "[directional_reduce]")
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // mesh
    int ngpts[3] = { 128, 256, 196 };

    // mesh spacing
    double hh = 0.1;

    // distance for data gathering
    double radius = 4.;

    pb::PEenv myPEenv(MPI_COMM_WORLD, ngpts[0], ngpts[1], ngpts[2], 1);

    double domain[3] = { ngpts[0] * hh, ngpts[1] * hh, ngpts[2] * hh };

    DirectionalReduce dir_reduce(myPEenv.cart_comm(), radius, domain);

    int nprocs[3];
    int periodic[3];
    int coords[3];
    // int error =
    MPI_Cart_get(myPEenv.cart_comm(), 3, nprocs, periodic, coords);

    // direction for test
    const short dir = 1;

    if (myrank == 0)
        std::cout << "nprocs[" << dir << "]=" << nprocs[dir] << std::endl;
    if (myrank == 0)
        std::cout << "rstep[" << dir << "]=" << dir_reduce.rstep(dir)
                  << std::endl;

    // compute max. coord in direction dir with radius 'radius'
    int data[2] = { coords[dir], 1 };
    dir_reduce.computeDirMax(dir, data);

    int result = std::min(coords[dir] + dir_reduce.rstep(dir), nprocs[dir] - 1);
    if (dir_reduce.lstep(dir) > coords[dir]) result = nprocs[dir] - 1;

    std::cout << "myrank = " << myrank << ", coords[dir] = " << coords[dir]
              << ", data[0] = " << data[0] << std::endl;

    CHECK(result == data[0]);
}
