// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "mgmol_mpi_tools.h"

#include <cassert>
#include <iostream>
#include <string.h>
using namespace std;

namespace mgmol_tools
{
int reduce(int* sendbuf, int* recvbuf, int count, MPI_Op op, const int root,
    const MPI_Comm comm)
{
#ifdef USE_MPI
    assert(comm != 0);
    int mpi_err = MPI_Reduce(sendbuf, recvbuf, count, MPI_INT, op, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Reduce(int*, int*) of size " << count << "!!!"
             << endl;
    }
    return mpi_err;
#else
    return 0;
#endif
}

int reduce(int* buf, int count, MPI_Op op, const int root, const MPI_Comm comm)
{
    int mype = 0;
#ifdef USE_MPI
    assert(comm != 0);
    MPI_Comm_rank(comm, &mype);
#endif

    int* recvbuf = new int[count];
    int ret      = reduce(buf, recvbuf, count, op, root, comm);

    if (mype == root) memcpy(buf, recvbuf, count * sizeof(int));
    delete[] recvbuf;

    return ret;
}

int reduce(double* sendbuf, double* recvbuf, int count, MPI_Op op,
    const int root, const MPI_Comm comm)
{
#ifdef USE_MPI
    assert(comm != 0);
    int mpi_err
        = MPI_Reduce(sendbuf, recvbuf, count, MPI_DOUBLE, op, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Reduce(double*, double*) of size " << count
             << "!!!" << endl;
    }
    return mpi_err;
#else
    return 0;
#endif
}

int reduce(float* sendbuf, float* recvbuf, int count, MPI_Op op, const int root,
    const MPI_Comm comm)
{
#ifdef USE_MPI
    assert(comm != 0);
    int mpi_err
        = MPI_Reduce(sendbuf, recvbuf, count, MPI_FLOAT, op, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Reduce(float*, float*) of size " << count << "!!!"
             << endl;
    }
    return mpi_err;
#else
    return 0;
#endif
}

int reduce(short* sendbuf, short* recvbuf, int count, MPI_Op op, const int root,
    const MPI_Comm comm)
{
#ifdef USE_MPI
    assert(comm != 0);
    int mpi_err
        = MPI_Reduce(sendbuf, recvbuf, count, MPI_SHORT, op, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Reduce(short*, short*) of size " << count << "!!!"
             << endl;
    }
    return mpi_err;
#else
    return 0;
#endif
}

int gatherV(vector<string>& sendbuf, vector<string>& recvbuf, const int root,
    const MPI_Comm comm)
{
    int mype = 0;
    int size = 1;
#ifdef USE_MPI
    assert(comm != 0);
    MPI_Comm_rank(comm, &mype);
    MPI_Comm_size(comm, &size);

    int vcount = (int)sendbuf.size();
    reduce(&vcount, 1, MPI_SUM, root, comm);

    int sendcount;
    int totchars = 0;
    // first get length of each string
    vector<int> locstrlen;
    vector<int> strLen(vcount, 0);
    sendcount = 0;
    for (vector<string>::iterator str = sendbuf.begin(); str != sendbuf.end();
         ++str)
    {
        string s = *str;
        sendcount += s.length();
        locstrlen.push_back(s.length());
    }
    gatherV(locstrlen, strLen, root, comm);

    // convert string to char array
    //    vector<char>charStr(sendcount);
    char* charStr = new char[sendcount];
    int idx       = 0;
    for (vector<string>::iterator str = sendbuf.begin(); str != sendbuf.end();
         ++str)
    {
        string s = *str;
        memcpy(&charStr[idx], s.c_str(), s.size());
        idx += s.size();
    }
    assert(idx == sendcount);

    int* recvcounts = new int[size];
    int mpi_err     = MPI_Gather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, root, comm);

    int* displs = new int[size];
    if (mype == root)
    {
        displs[0] = 0;
        totchars  = recvcounts[size - 1];
        for (int i = 1; i < size; i++)
        {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
            totchars += recvcounts[i - 1];
        }
    }
    char* recvdata = new char[totchars];
    mpi_err        = MPI_Gatherv(&charStr[0], sendcount, MPI_CHAR, &recvdata[0],
        recvcounts, displs, MPI_CHAR, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Gatherv in MGmol_MPI::GatherV() !!!" << endl;
    }

    // reset recvbuf
    recvbuf.clear();
    if (mype == root)
    {
        int pos = 0;
        for (int i = 0; i < vcount; i++)
        {
            char str[strLen[i] + 1];
            str[strLen[i]] = '\0';
            memcpy(str, &recvdata[pos], strLen[i] * sizeof(char));
            string cstr;
            cstr.assign(str);
            recvbuf.push_back(cstr);
            pos += strLen[i] * sizeof(char);
        }
    }

    delete[] displs;
    delete[] recvcounts;
    delete[] charStr;
    delete[] recvdata;

    return mpi_err;
#else
    return 0;
#endif
}

int gatherV(vector<int>& sendbuf, vector<int>& recvbuf, const int root,
    const MPI_Comm comm)
{
#ifdef USE_MPI
    int mype = 0;
    int size = 1;
    MPI_Comm_rank(comm, &mype);
    MPI_Comm_size(comm, &size);

    int sendcount   = (int)sendbuf.size();
    int* recvcounts = new int[size];
    int mpi_err     = MPI_Gather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Allgather in MGmol_MPI::gatherV() !!!" << endl;
    }

    vector<int> displs(size);
    if (mype == root)
    {
        displs[0] = 0;
        for (int i = 1; i < size; i++)
        {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
        }

        recvbuf.resize(displs[size - 1] + recvcounts[size - 1]);
    }
    mpi_err = MPI_Gatherv(&sendbuf[0], sendcount, MPI_INT, &recvbuf[0],
        recvcounts, &displs[0], MPI_INT, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Gatherv in MGmol_MPI::gatherV() !!!" << endl;
    }

    delete[] recvcounts;

    return mpi_err;
#else
    return 0;
#endif
}

int gatherV(vector<unsigned short>& sendbuf, vector<unsigned short>& recvbuf,
    const int root, const MPI_Comm comm)
{
#ifdef USE_MPI
    int mype = 0;
    int size = 1;
    MPI_Comm_rank(comm, &mype);
    MPI_Comm_size(comm, &size);

    int sendcount   = (int)sendbuf.size();
    int* recvcounts = new int[size];
    int mpi_err     = MPI_Gather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Allgather in MGmol_MPI::gatherV() !!!" << endl;
    }

    int* displs = new int[size];
    if (mype == root)
    {
        displs[0] = 0;
        for (int i = 1; i < size; i++)
        {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
        }

        recvbuf.resize(displs[size - 1] + recvcounts[size - 1]);
    }
    mpi_err = MPI_Gatherv(&sendbuf[0], sendcount, MPI_UNSIGNED_SHORT,
        &recvbuf[0], recvcounts, displs, MPI_UNSIGNED_SHORT, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Gatherv in MGmol_MPI::gatherV() !!!" << endl;
    }

    delete[] displs;
    delete[] recvcounts;

    return mpi_err;
#else
    return 0;
#endif
}

int gatherV(vector<double>& sendbuf, vector<double>& recvbuf, const int root,
    const MPI_Comm comm)
{
#ifdef USE_MPI
    int mype = 0;
    int size = 1;
    MPI_Comm_rank(comm, &mype);
    MPI_Comm_size(comm, &size);

    int sendcount   = (int)sendbuf.size();
    int* recvcounts = new int[size];
    int mpi_err     = MPI_Gather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Allgather in MGmol_MPI::gatherV() !!!" << endl;
    }

    int* displs = new int[size];
    if (mype == root)
    {
        displs[0] = 0;
        for (int i = 1; i < size; i++)
        {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
        }

        recvbuf.resize(displs[size - 1] + recvcounts[size - 1]);
    }
    mpi_err = MPI_Gatherv(&sendbuf[0], sendcount, MPI_DOUBLE, &recvbuf[0],
        recvcounts, displs, MPI_DOUBLE, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Gatherv in MGmol_MPI::gatherV() !!!" << endl;
    }

    delete[] displs;
    delete[] recvcounts;

    return mpi_err;
#else
    return 0;
#endif
}

int gatherV(vector<float>& sendbuf, vector<float>& recvbuf, const int root,
    const MPI_Comm comm)
{
#ifdef USE_MPI
    int mype = 0;
    int size = 1;
    MPI_Comm_rank(comm, &mype);
    MPI_Comm_size(comm, &size);

    int sendcount   = (int)sendbuf.size();
    int* recvcounts = new int[size];
    int mpi_err     = MPI_Gather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Allgather in MGmol_MPI::gatherV() !!!" << endl;
    }

    int* displs = new int[size];
    if (mype == root)
    {
        displs[0] = 0;
        for (int i = 1; i < size; i++)
        {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
        }

        recvbuf.resize(displs[size - 1] + recvcounts[size - 1]);
    }
    mpi_err = MPI_Gatherv(&sendbuf[0], sendcount, MPI_FLOAT, &recvbuf[0],
        recvcounts, displs, MPI_FLOAT, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Gatherv in MGmol_MPI::gatherV() !!!" << endl;
    }

    delete[] displs;
    delete[] recvcounts;

    return mpi_err;
#else
    return 0;
#endif
}

int allreduce(
    short* sendbuf, short* recvbuf, int count, MPI_Op op, const MPI_Comm comm)
{
#ifdef USE_MPI
    int mpi_err = MPI_Allreduce(sendbuf, recvbuf, count, MPI_SHORT, op, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        cerr << "ERROR in MPI_Allreduce(int*, int*) of size " << count << "!!!"
             << endl;
    }
    return mpi_err;
#else
    return 0;
#endif
}
};
