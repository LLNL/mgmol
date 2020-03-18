// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "mpi.h"
#include <cassert>
#include <iostream>
#include <string.h>

void MPI_Irecv(double* precv, const int blocksize, const int src, const int i,
    MPI_Comm comm, MPI_Request& rr)
{
    MPI_Irecv(precv, blocksize, MPI_DOUBLE, src, i, comm, &rr);
}

void MPI_Irecv(int* precv, const int blocksize, const int src, const int i,
    MPI_Comm comm, MPI_Request& rr)
{
    MPI_Irecv(precv, blocksize, MPI_INT, src, i, comm, &rr);
}

void MPI_Irecv(unsigned short* precv, const int blocksize, const int src,
    const int i, MPI_Comm comm, MPI_Request& rr)
{
    MPI_Irecv(precv, blocksize, MPI_UNSIGNED_SHORT, src, i, comm, &rr);
}

void MPI_Isend(double* psend, const int blocksize, const int dst, const int i,
    MPI_Comm comm)
{
    MPI_Request sr;
    MPI_Isend(psend, blocksize, MPI_DOUBLE, dst, i, comm, &sr);
}

void MPI_Isend(
    int* psend, const int blocksize, const int dst, const int i, MPI_Comm comm)
{
    MPI_Request sr;
    MPI_Isend(psend, blocksize, MPI_INT, dst, i, comm, &sr);
}

void MPI_Isend(unsigned short* psend, const int blocksize, const int dst,
    const int i, MPI_Comm comm)
{
    MPI_Request sr;
    MPI_Isend(psend, blocksize, MPI_UNSIGNED_SHORT, dst, i, comm, &sr);
}

template <typename ScalarType>
int MPI_AlltofirstN(ScalarType* sendbuf, const int blocksize,
    ScalarType* recvbuf, const int nfirst, MPI_Comm comm)
{

    int npes_send;
    MPI_Comm_size(comm, &npes_send);

    int src;
    MPI_Comm_rank(comm, &src);
    int dst         = src;
    const bool recv = (src < nfirst);
    if (recv) assert(recvbuf != nullptr);

    if (recv)
        memcpy(recvbuf + src * blocksize, sendbuf + dst * blocksize,
            blocksize * sizeof(ScalarType));

    for (int i = 1; i < npes_send; i++)
    {
        MPI_Request rr;

        // increment source
        src++;
        if (src >= npes_send) src = 0;

        dst--;
        if (dst < 0) dst = npes_send - 1;

        // recv data from src
        if (recv)
        {
            ScalarType* precv = recvbuf + src * blocksize;
            MPI_Irecv(precv, blocksize, src, i, comm, rr);
        }
        if (dst < nfirst)
        {
            ScalarType* psend = sendbuf + dst * blocksize;
            MPI_Isend(psend, blocksize, dst, i, comm);
        }

        if (recv) MPI_Wait(&rr, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(comm);

    return 0;
}

template <typename ScalarType>
int MPI_AlltofirstNv(ScalarType* sendbuf, int* sendcnts, int* sdispls,
    ScalarType* recvbuf, int* recvcnts, int* rdispls, const int nfirst,
    MPI_Comm comm)
{
    MPI_Request rr;

    int npes_send;
    MPI_Comm_size(comm, &npes_send);

    int src;
    MPI_Comm_rank(comm, &src);
    int dst         = src;
    const bool recv = (src < nfirst);
    if (recv) assert(recvbuf != nullptr);

    if (recv)
        memcpy(recvbuf + rdispls[src], sendbuf + sdispls[dst],
            recvcnts[src] * sizeof(ScalarType));

    for (int i = 1; i < npes_send; i++)
    {
        // increment source
        src++;
        if (src >= npes_send) src = 0;

        dst--;
        if (dst < 0) dst = npes_send - 1;

        // recv data from src
        if (recv)
        {
            ScalarType* precv = recvbuf + rdispls[src];
            if (recvcnts[src] > 0)
                MPI_Irecv(precv, recvcnts[src], src, i, comm, rr);
        }
        if (dst < nfirst)
        {
            ScalarType* psend = sendbuf + sdispls[dst];
            if (sendcnts[dst] > 0)
                MPI_Isend(psend, sendcnts[dst], dst, i, comm);
        }
        if (recv && recvcnts[src] > 0) MPI_Wait(&rr, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(comm);

    return 0;
}

// Explicit instantiation
template int MPI_AlltofirstN<double>(double* sendbuf, const int blocksize,
    double* recvbuf, const int nfirst, MPI_Comm comm);

template int MPI_AlltofirstN<int>(int* sendbuf, const int blocksize,
    int* recvbuf, const int nfirst, MPI_Comm comm);

template int MPI_AlltofirstN<unsigned short>(unsigned short* sendbuf,
    const int blocksize, unsigned short* recvbuf, const int nfirst,
    MPI_Comm comm);

template int MPI_AlltofirstNv<int>(int* sendbuf, int* sendcnts, int* sdispls,
    int* recvbuf, int* recvcnts, int* rdispls, const int nfirst, MPI_Comm comm);

template int MPI_AlltofirstNv<unsigned short>(unsigned short* sendbuf,
    int* sendcnts, int* sdispls, unsigned short* recvbuf, int* recvcnts,
    int* rdispls, const int nfirst, MPI_Comm comm);
