// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_MPI_TOOLS_H
#define MGMOL_MPI_TOOLS_H

#include <mpi.h>
#include <string>
#include <vector>

namespace mgmol_tools
{
int reduce(int* sendbuf, int* recvbuf, int count, MPI_Op op, const int root,
    const MPI_Comm comm);
int reduce(double* sendbuf, double* recvbuf, int count, MPI_Op op,
    const int root, const MPI_Comm comm);
int reduce(float* sendbuf, float* recvbuf, int count, MPI_Op op, const int root,
    const MPI_Comm comm);
int reduce(short* sendbuf, short* recvbuf, int count, MPI_Op op, const int root,
    const MPI_Comm comm);

int reduce(int* buf, int count, MPI_Op op, const int root, const MPI_Comm comm);

int gatherV(std::vector<std::string>& sendbuf,
    std::vector<std::string>& recvbuf, const int root, const MPI_Comm comm);
int gatherV(std::vector<int>& sendbuf, std::vector<int>& recvbuf,
    const int root, const MPI_Comm comm);
int gatherV(std::vector<double>& sendbuf, std::vector<double>& recvbuf,
    const int root, const MPI_Comm comm);
int gatherV(std::vector<float>& sendbuf, std::vector<float>& recvbuf,
    const int root, const MPI_Comm comm);
int gatherV(std::vector<unsigned short>& sendbuf,
    std::vector<unsigned short>& recvbuf, const int root, const MPI_Comm comm);

int allreduce(
    short* sendbuf, short* recvbuf, int count, MPI_Op op, const MPI_Comm comm);
}

#endif
