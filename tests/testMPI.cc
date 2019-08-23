// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include <iostream>
#include <mpi.h>

int main(int argc, char** argv)
{
    /*
     * Initialize MPI
     */
    int mpi_ret = MPI_Init(&argc, &argv);
    if (mpi_ret != MPI_SUCCESS)
    {
        std::cerr << "Error in MPI_Init..." << std::endl;
        return 1;
    }

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs != 4)
    {
        std::cerr << "Expected 4 MPI ranks. Got " << nprocs << std::endl;
        return 1;
    }

    int myproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

    int sendbuff = myproc;
    int recvbuf  = 0;
    int sum      = nprocs * (nprocs + 1) / 2 - nprocs;
    mpi_ret      = MPI_Allreduce(
        &sendbuff, &recvbuf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (mpi_ret == MPI_ERR_COMM)
    {
        std::cerr << "Invalid communicator!" << std::endl;
        return 1;
    }
    if (recvbuf != sum)
    {
        std::cerr << "Incorrect sum for myproc()!!" << std::endl;
        return 1;
    }

    mpi_ret = MPI_Finalize();
    if (mpi_ret != MPI_SUCCESS)
    {
        std::cerr << "Error in MPI_Finalize..." << std::endl;
        ;
        return 1;
    }

    std::cout << "SUCESSFUL test!" << std::endl;

    return 0;
}
