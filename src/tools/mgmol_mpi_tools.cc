// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "mgmol_mpi_tools.h"

#include <cassert>
#include <iostream>
#include <memory>
#include <string.h>

namespace mgmol_tools
{
int reduce(int* sendbuf, int* recvbuf, int count, MPI_Op op, const int root,
    const MPI_Comm comm)
{
    assert(comm != MPI_COMM_NULL);
    int mpi_err = MPI_Reduce(sendbuf, recvbuf, count, MPI_INT, op, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        std::cerr << "ERROR in MPI_Reduce(int*, int*) of size " << count
                  << "!!!" << std::endl;
    }
    return mpi_err;
}

int reduce(int* buf, int count, MPI_Op op, const int root, const MPI_Comm comm)
{
    int mype = 0;
    assert(comm != MPI_COMM_NULL);
    MPI_Comm_rank(comm, &mype);

    int* recvbuf = new int[count];
    int ret      = reduce(buf, recvbuf, count, op, root, comm);

    if (mype == root) memcpy(buf, recvbuf, count * sizeof(int));
    delete[] recvbuf;

    return ret;
}

int reduce(double* sendbuf, double* recvbuf, int count, MPI_Op op,
    const int root, const MPI_Comm comm)
{
    assert(comm != MPI_COMM_NULL);
    int mpi_err
        = MPI_Reduce(sendbuf, recvbuf, count, MPI_DOUBLE, op, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        std::cerr << "ERROR in MPI_Reduce(double*, double*) of size " << count
                  << "!!!" << std::endl;
    }
    return mpi_err;
}

int reduce(float* sendbuf, float* recvbuf, int count, MPI_Op op, const int root,
    const MPI_Comm comm)
{
    assert(comm != MPI_COMM_NULL);
    int mpi_err
        = MPI_Reduce(sendbuf, recvbuf, count, MPI_FLOAT, op, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        std::cerr << "ERROR in MPI_Reduce(float*, float*) of size " << count
                  << "!!!" << std::endl;
    }
    return mpi_err;
}

int reduce(short* sendbuf, short* recvbuf, int count, MPI_Op op, const int root,
    const MPI_Comm comm)
{
    assert(comm != MPI_COMM_NULL);
    int mpi_err
        = MPI_Reduce(sendbuf, recvbuf, count, MPI_SHORT, op, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        std::cerr << "ERROR in MPI_Reduce(short*, short*) of size " << count
                  << "!!!" << std::endl;
    }
    return mpi_err;
}

int gatherV(std::vector<std::string>& sendbuf,
    std::vector<std::string>& recvbuf, const int root, const MPI_Comm comm)
{
    int mype = 0;
    int size = 1;
    assert(comm != MPI_COMM_NULL);
    MPI_Comm_rank(comm, &mype);
    MPI_Comm_size(comm, &size);

    int vcount = (int)sendbuf.size();
    reduce(&vcount, 1, MPI_SUM, root, comm);

    int sendcount;
    int totchars = 0;
    // first get length of each string
    std::vector<int> locstrlen;
    std::vector<int> strLen(vcount, 0);
    sendcount = 0;
    for (std::vector<std::string>::iterator str = sendbuf.begin();
         str != sendbuf.end(); ++str)
    {
        std::string s = *str;
        sendcount += s.length();
        locstrlen.push_back(s.length());
    }
    gatherV(locstrlen, strLen, root, comm);

    // convert string to char array
    char* charStr = new char[sendcount];
    int idx       = 0;
    for (std::vector<std::string>::iterator str = sendbuf.begin();
         str != sendbuf.end(); ++str)
    {
        std::string s = *str;
        memcpy(&charStr[idx], s.c_str(), s.size());
        idx += s.size();
    }
    assert(idx == sendcount);

    int* recvcounts = new int[size];
    int mpi_err     = MPI_Gather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, root, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        std::cerr << "ERROR in MPI_Gather in MGmol_MPI::GatherV() !!!"
                  << std::endl;
    }

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
        std::cerr << "ERROR in MPI_Gatherv in MGmol_MPI::GatherV() !!!"
                  << std::endl;
    }

    // reset recvbuf
    recvbuf.clear();
    if (mype == root)
    {
        int pos = 0;
        for (int i = 0; i < vcount; i++)
        {
            std::unique_ptr<char[]> str(new char[strLen[i] + 1]);
            str[strLen[i]] = '\0';
            memcpy(str.get(), &recvdata[pos], strLen[i] * sizeof(char));
            std::string cstr;
            cstr.assign(str.get());
            recvbuf.push_back(cstr);
            pos += strLen[i] * sizeof(char);
        }
    }

    delete[] displs;
    delete[] recvcounts;
    delete[] charStr;
    delete[] recvdata;

    return mpi_err;
}

int gatherV(std::vector<int>& sendbuf, std::vector<int>& recvbuf,
    const int root, const MPI_Comm comm)
{
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
        std::cerr << "ERROR in MPI_Allgather in MGmol_MPI::gatherV() !!!"
                  << std::endl;
    }

    std::vector<int> displs(size);
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
        std::cerr << "ERROR in MPI_Gatherv in MGmol_MPI::gatherV() !!!"
                  << std::endl;
    }

    delete[] recvcounts;

    return mpi_err;
}

int gatherV(std::vector<unsigned short>& sendbuf,
    std::vector<unsigned short>& recvbuf, const int root, const MPI_Comm comm)
{
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
        std::cerr << "ERROR in MPI_Allgather in MGmol_MPI::gatherV() !!!"
                  << std::endl;
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
        std::cerr << "ERROR in MPI_Gatherv in MGmol_MPI::gatherV() !!!"
                  << std::endl;
    }

    delete[] displs;
    delete[] recvcounts;

    return mpi_err;
}

int gatherV(std::vector<double>& sendbuf, std::vector<double>& recvbuf,
    const int root, const MPI_Comm comm)
{
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
        std::cerr << "ERROR in MPI_Allgather in MGmol_MPI::gatherV() !!!"
                  << std::endl;
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
        std::cerr << "ERROR in MPI_Gatherv in MGmol_MPI::gatherV() !!!"
                  << std::endl;
    }

    delete[] displs;
    delete[] recvcounts;

    return mpi_err;
}

int gatherV(std::vector<float>& sendbuf, std::vector<float>& recvbuf,
    const int root, const MPI_Comm comm)
{
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
        std::cerr << "ERROR in MPI_Allgather in MGmol_MPI::gatherV() !!!"
                  << std::endl;
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
        std::cerr << "ERROR in MPI_Gatherv in MGmol_MPI::gatherV() !!!"
                  << std::endl;
    }

    delete[] displs;
    delete[] recvcounts;

    return mpi_err;
}

int allreduce(
    short* sendbuf, short* recvbuf, int count, MPI_Op op, const MPI_Comm comm)
{
    int mpi_err = MPI_Allreduce(sendbuf, recvbuf, count, MPI_SHORT, op, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        std::cerr << "ERROR in MPI_Allreduce(int*, int*) of size " << count
                  << "!!!" << std::endl;
    }
    return mpi_err;
}
}
