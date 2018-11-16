// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: MPI_AllToSome.C,v 1.3 2010/09/29 22:48:08 jeanluc Exp $
#ifdef USE_MPI

#include "mpi.h"
#include <cassert>
#include <iostream>
#include <string.h>
using namespace std;

int MPI_AlltofirstN(double* sendbuf, const int blocksize, double* recvbuf,
    const int nfirst, MPI_Comm comm)
{
    MPI_Request rr;
    MPI_Request sr;

    int npes_send;
    MPI_Comm_size(comm, &npes_send);

    int src;
    MPI_Comm_rank(comm, &src);
    int dst         = src;
    const bool recv = (src < nfirst);
    if (recv) assert(recvbuf != NULL);

    if (recv)
        memcpy(recvbuf + src * blocksize, sendbuf + dst * blocksize,
            blocksize * sizeof(double));

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
            double* precv = recvbuf + src * blocksize;
            MPI_Irecv(precv, blocksize, MPI_DOUBLE, src, i, comm, &rr);
        }
        if (dst < nfirst)
        {
            double* psend = sendbuf + dst * blocksize;
            MPI_Isend(psend, blocksize, MPI_DOUBLE, dst, i, comm, &sr);
        }

        if (recv) MPI_Wait(&rr, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(comm);

    return 0;
}

int MPI_AlltofirstN(int* sendbuf, const int blocksize, int* recvbuf,
    const int nfirst, MPI_Comm comm)
{
    MPI_Request rr;
    MPI_Request sr;

    int npes_send;
    MPI_Comm_size(comm, &npes_send);

    int src;
    MPI_Comm_rank(comm, &src);
    int dst         = src;
    const bool recv = (src < nfirst);
    if (recv) assert(recvbuf != NULL);
    // if( src==0 )cout<<"MPI_AlltofirstN with nfirst="<<nfirst<<endl;

    if (recv)
        memcpy(recvbuf + src * blocksize, sendbuf + dst * blocksize,
            blocksize * sizeof(int));

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
            int* precv = recvbuf + src * blocksize;
            MPI_Irecv(precv, blocksize, MPI_INT, src, i, comm, &rr);
        }
        if (dst < nfirst)
        {
            int* psend = sendbuf + dst * blocksize;
            MPI_Isend(psend, blocksize, MPI_INT, dst, i, comm, &sr);
        }

        if (recv) MPI_Wait(&rr, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(comm);

    return 0;
}

int MPI_AlltofirstN(unsigned short* sendbuf, const int blocksize,
    unsigned short* recvbuf, const int nfirst, MPI_Comm comm)
{
    MPI_Request rr;
    MPI_Request sr;

    int npes_send;
    MPI_Comm_size(comm, &npes_send);

    int src;
    MPI_Comm_rank(comm, &src);
    int dst         = src;
    const bool recv = (src < nfirst);
    if (recv) assert(recvbuf != NULL);
    // if( src==0 )cout<<"MPI_AlltofirstN with nfirst="<<nfirst<<endl;

    if (recv)
        memcpy(recvbuf + src * blocksize, sendbuf + dst * blocksize,
            blocksize * sizeof(unsigned short));

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
            unsigned short* precv = recvbuf + src * blocksize;
            MPI_Irecv(precv, blocksize, MPI_UNSIGNED_SHORT, src, i, comm, &rr);
        }
        if (dst < nfirst)
        {
            unsigned short* psend = sendbuf + dst * blocksize;
            MPI_Isend(psend, blocksize, MPI_UNSIGNED_SHORT, dst, i, comm, &sr);
        }

        if (recv) MPI_Wait(&rr, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(comm);

    return 0;
}

int MPI_AlltofirstNv(int* sendbuf, int* sendcnts, int* sdispls, int* recvbuf,
    int* recvcnts, int* rdispls, const int nfirst, MPI_Comm comm)
{
    MPI_Request rr;
    MPI_Request sr;

    int npes_send;
    MPI_Comm_size(comm, &npes_send);

    int src;
    MPI_Comm_rank(comm, &src);
    int dst         = src;
    const bool recv = (src < nfirst);
    if (recv) assert(recvbuf != NULL);
    // if( src==0 )cout<<"MPI_AlltofirstNv with nfirst="<<nfirst<<endl;

    if (recv)
        memcpy(recvbuf + rdispls[src], sendbuf + sdispls[dst],
            recvcnts[src] * sizeof(int));

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
            int* precv = recvbuf + rdispls[src];
            if (recvcnts[src] > 0)
                MPI_Irecv(precv, recvcnts[src], MPI_INT, src, i, comm, &rr);
        }
        if (dst < nfirst)
        {
            int* psend = sendbuf + sdispls[dst];
            if (sendcnts[dst] > 0)
                MPI_Isend(psend, sendcnts[dst], MPI_INT, dst, i, comm, &sr);
        }
        if (recv && recvcnts[src] > 0) MPI_Wait(&rr, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(comm);

    return 0;
}

int MPI_AlltofirstNv(unsigned short* sendbuf, int* sendcnts, int* sdispls,
    unsigned short* recvbuf, int* recvcnts, int* rdispls, const int nfirst,
    MPI_Comm comm)
{
    MPI_Request rr;
    MPI_Request sr;

    int npes_send;
    MPI_Comm_size(comm, &npes_send);

    int src;
    MPI_Comm_rank(comm, &src);
    int dst         = src;
    const bool recv = (src < nfirst);
    if (recv) assert(recvbuf != NULL);
    // if( src==0 )cout<<"MPI_AlltofirstNv with nfirst="<<nfirst<<endl;

    if (recv)
        memcpy(recvbuf + rdispls[src], sendbuf + sdispls[dst],
            recvcnts[src] * sizeof(unsigned short));

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
            unsigned short* precv = recvbuf + rdispls[src];
            if (recvcnts[src] > 0)
                MPI_Irecv(precv, recvcnts[src], MPI_UNSIGNED_SHORT, src, i,
                    comm, &rr);
        }
        if (dst < nfirst)
        {
            unsigned short* psend = sendbuf + sdispls[dst];
            if (sendcnts[dst] > 0)
                MPI_Isend(psend, sendcnts[dst], MPI_UNSIGNED_SHORT, dst, i,
                    comm, &sr);
        }
        if (recv && recvcnts[src] > 0) MPI_Wait(&rr, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(comm);

    return 0;
}

#endif
