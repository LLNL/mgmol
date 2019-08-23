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

#include <iostream>

int main(int argc, char** argv)
{
    /*
     * Initialize MPI
     */
    int mpi_ret = MPI_Init(&argc, &argv);
    if (mpi_ret != MPI_SUCCESS)
    {
        std::cerr << "testBlacsContext: Error in MPI_Init..." << std::endl;
        return 1;
    }

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    {

        // build a BlacsContext 1 row x nprocs col
        dist_matrix::BlacsContext bc(MPI_COMM_WORLD, 'r', nprocs);

        if (bc.nprow() != 1)
        {
            std::cerr << "Wrong number of rows!!" << std::endl;
            return 1;
        }
        if (bc.npcol() != nprocs)
        {
            std::cerr << "Wrong number of columns!!" << std::endl;
            return 1;
        }
        if (bc.nprocs() != nprocs)
        {
            std::cerr << "Wrong number of processors!!" << std::endl;
            return 1;
        }

        if (bc.myproc() >= nprocs || bc.myproc() < 0)
        {
            std::cerr << "Incorrect value for myproc()!!" << std::endl;
            return 1;
        }

        int sendbuff = bc.myproc();
        int recvbuf  = 0;
        int sum      = nprocs * (nprocs + 1) / 2 - nprocs;
        MPI_Allreduce(&sendbuff, &recvbuf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (recvbuf != sum)
        {
            std::cerr << "Incorrect sum for myproc()!!" << std::endl;
            return 1;
        }
    }

    mpi_ret = MPI_Finalize();
    if (mpi_ret != MPI_SUCCESS)
    {
        std::cerr << "testBlacsContext: Error in MPI_Finalize..." << std::endl;
        ;
        return 1;
    }

    return 0;
}
