// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "Timer.h"
#include "mgmol_mpi_tools.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string.h>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX_SIZE (262144)

Timer MGmol_MPI::split_allreduce_sums_double_tm_("split_allreduce_sums_double");
Timer MGmol_MPI::split_allreduce_sums_float_tm_("split_allreduce_sums_float");

MGmol_MPI* MGmol_MPI::pinstance_ = nullptr;
std::ostream* MGmol_MPI::os_     = nullptr;

MPI_Comm MGmol_MPI::comm_global_         = MPI_COMM_NULL;
MPI_Comm MGmol_MPI::comm_spin_           = MPI_COMM_NULL;
MPI_Comm MGmol_MPI::comm_different_spin_ = MPI_COMM_NULL;
MPI_Comm MGmol_MPI::comm_images_         = MPI_COMM_NULL;

int MGmol_MPI::mype_      = -1;
int MGmol_MPI::mype_spin_ = -1;
int MGmol_MPI::size_      = -1;

short MGmol_MPI::myspin_ = -1;
short MGmol_MPI::nspin_  = -1;

short MGmol_MPI::myimage_ = -1;
short MGmol_MPI::nimages_ = -1;

MGmol_MPI::MGmol_MPI()
{
    assert(comm_spin_ != MPI_COMM_NULL);

    MPI_Comm_size(comm_spin_, &size_);
    MPI_Comm_rank(comm_spin_, &mype_spin_);
    if (mype_spin_ == 0)
        *os_ << "MGmol instance using " << size_ << " MPI tasks" << std::endl;
}

void MGmol_MPI::printTimers(std::ostream& os)
{
    split_allreduce_sums_double_tm_.print(os);
    split_allreduce_sums_float_tm_.print(os);
}

void MGmol_MPI::setupComm(
    const MPI_Comm comm, const bool with_spin, const int nimages)
{
    assert(pinstance_ == nullptr);

    comm_global_ = comm;

    int npes;
    MPI_Comm_rank(comm, &mype_);
    MPI_Comm_size(comm, &npes);

    nimages_ = nimages;

    if (nimages_ > 1)
    {
        // cout<<"nimages="<<nimages_<<endl;
        // create communicator to communicate data within one image calculation
        int npes_sub = npes / nimages_;
        assert(nimages_ * npes_sub == npes);
        myimage_ = mype_ / npes_sub; // 0 to nimages_-1
        assert(myimage_ >= 0 && myimage_ < nimages_);
        int key = mype_ % npes_sub;
#ifndef NDEBUG
        int mpirc =
#endif
            MPI_Comm_split(comm, myimage_, key, &comm_spin_);
        // cout<<"myimage_="<<myimage_<<", key="<<key<<", comm_="<<comm_<<endl;
        assert(mpirc == MPI_SUCCESS);

        // create communicator to communicate data from other images
        int color = mype_ % npes_sub;
#ifndef NDEBUG
        mpirc =
#endif
            MPI_Comm_split(comm, color, myimage_, &comm_images_);
        assert(mpirc == MPI_SUCCESS);

        MPI_Comm_size(comm_images_, &npes);
        assert(npes == nimages_);

        nspin_  = 1;
        myspin_ = 0;
    }
    else
    {
        if (with_spin)
        {
            if (npes % 2 != 0 && mype_ == 0)
            {
                std::cerr
                    << " Calculation with spin requires even number of MPI "
                       "tasks!!!"
                    << std::endl;
                MPI_Abort(comm, EXIT_FAILURE);
            }

            // create communicator to communicate data within one spin
            // calculation
            int npes_sub = npes / 2;
            int color    = mype_ / npes_sub; // 0 or 1
            assert(color == 0 || color == 1);
            int key = mype_ % npes_sub;
            assert(key < npes_sub);
            int mpierr = MPI_Comm_split(comm_global_, color, key, &comm_spin_);
            if (mpierr != MPI_SUCCESS)
            {
                std::cerr << " Error in creating spin subcommunicator!!!"
                          << std::endl;
                MPI_Abort(comm, EXIT_FAILURE);
            }

            MPI_Barrier(comm_global_);

            // create communicator to communicate data from other spin
            color = mype_ % npes_sub;
            key   = mype_ / npes_sub; // 0 or 1
            assert(key == 0 || key == 1);
            assert(color < npes_sub);
            mpierr = MPI_Comm_split(
                comm_global_, color, key, &comm_different_spin_);
            if (mpierr != MPI_SUCCESS)
            {
                std::cerr << " Error in creating across spin subcommunicator!!!"
                          << std::endl;
                MPI_Abort(comm, EXIT_FAILURE);
            }
            nspin_  = 2;
            myspin_ = key;

            MPI_Comm_size(comm_different_spin_, &npes);
            assert(npes == 2);
        }
        else
        {
            comm_spin_ = comm_global_;
            nspin_     = 1;
            myspin_    = 0;
        }
    }
    MPI_Comm_rank(comm_spin_, &mype_spin_);
}

int MGmol_MPI::bcast(double* val, int size, int root) const
{
    int mpi_err = MPI_Bcast(val, size, MPI_DOUBLE, root, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR(
            "ERROR in MPI_Bcast(double*, int) of size " << size << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::bcast(float* val, int size, int root) const
{
    int mpi_err = MPI_Bcast(val, size, MPI_FLOAT, root, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR(
            "ERROR in MPI_Bcast(float*, int) of size " << size << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::bcast(int* val, int size, int root) const
{
    int mpi_err = MPI_Bcast(val, size, MPI_INT, root, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Bcast(int*, int) of size " << size << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::bcast(short* val, int size, int root) const
{
    int mpi_err = MPI_Bcast(val, size, MPI_SHORT, root, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Bcast(short*, int) of size " << size << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::bcastGlobal(double* val, int size, int root) const
{
    int mpi_err = MPI_Bcast(val, size, MPI_DOUBLE, root, comm_global_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Bcast(double*, int) of size " << size << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::bcastGlobal(short* val, int size, int root) const
{
    int mpi_err = MPI_Bcast(val, size, MPI_SHORT, root, comm_global_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Bcast(short*, int) of size " << size << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::bcastGlobal(int* val, int size, int root) const
{
    int mpi_err = MPI_Bcast(val, size, MPI_INT, root, comm_global_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Bcast(short*, int) of size " << size << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::bcast(char* val, int size, int root) const
{
    int mpi_err = MPI_Bcast(val, size, MPI_CHAR, root, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Bcast(char*, int) of size " << size << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::reduce(double* sendbuf, double* recvbuf, int count, MPI_Op op,
    const int root) const
{
    return mgmol_tools::reduce(sendbuf, recvbuf, count, op, root, comm_spin_);
}

int MGmol_MPI::reduce(
    float* sendbuf, float* recvbuf, int count, MPI_Op op, const int root) const
{
    return mgmol_tools::reduce(sendbuf, recvbuf, count, op, root, comm_spin_);
}

int MGmol_MPI::reduce(
    int* sendbuf, int* recvbuf, int count, MPI_Op op, const int root) const
{
    return mgmol_tools::reduce(sendbuf, recvbuf, count, op, root, comm_spin_);
}

int MGmol_MPI::reduce(
    short* sendbuf, short* recvbuf, int count, MPI_Op op, const int root) const
{
    return mgmol_tools::reduce(sendbuf, recvbuf, count, op, root, comm_spin_);
}

int MGmol_MPI::allreduce(
    double* sendbuf, double* recvbuf, int count, MPI_Op op) const
{
    int mpi_err
        = MPI_Allreduce(sendbuf, recvbuf, count, MPI_DOUBLE, op, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR(
            "MPI_Allreduce(double*, double*) of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::allreduce(
    float* sendbuf, float* recvbuf, int count, MPI_Op op) const
{
    int mpi_err
        = MPI_Allreduce(sendbuf, recvbuf, count, MPI_FLOAT, op, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR(
            "MPI_Allreduce(float*, float*) of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::allreduce(int* sendbuf, int* recvbuf, int count, MPI_Op op) const
{
    int mpi_err
        = MPI_Allreduce(sendbuf, recvbuf, count, MPI_INT, op, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allreduce(int*, int*) of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::allreduce(
    short* sendbuf, short* recvbuf, int count, MPI_Op op) const
{
    int mpi_err
        = MPI_Allreduce(sendbuf, recvbuf, count, MPI_SHORT, op, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allreduce(int*, int*) of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::allreduceGlobal(
    int* sendbuf, int* recvbuf, int count, MPI_Op op) const
{
    int mpi_err
        = MPI_Allreduce(sendbuf, recvbuf, count, MPI_INT, op, comm_global_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allreduce(int*, int*) of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::allreduceGlobal(
    double* sendbuf, double* recvbuf, int count, MPI_Op op) const
{
    int mpi_err
        = MPI_Allreduce(sendbuf, recvbuf, count, MPI_DOUBLE, op, comm_global_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allreduce(int*, int*) of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::getRankMinVal(const double val, int* rank)
{
    typedef struct
    {
        double val;
        int rank;
    } VALRANK;

    VALRANK in;
    VALRANK out;
    in.val  = val;
    in.rank = mype_spin_;

    MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, comm_spin_);

    *rank       = out.rank;
    int mpi_err = MPI_Bcast(rank, 1, MPI_INT, 0, comm_spin_);

    return mpi_err;
}

int MGmol_MPI::getRankMaxVal(const float val, int* rank)
{
    typedef struct
    {
        float val;
        int rank;
    } VALRANK;

    VALRANK in;
    VALRANK out;
    in.val  = val;
    in.rank = mype_spin_;

    MPI_Reduce(&in, &out, 1, MPI_FLOAT_INT, MPI_MAXLOC, 0, comm_spin_);

    *rank       = out.rank;
    int mpi_err = MPI_Bcast(rank, 1, MPI_INT, 0, comm_spin_);

    return mpi_err;
}

int MGmol_MPI::allreduceImages(
    double* sendbuf, double* recvbuf, int count, MPI_Op op) const
{
    int mpi_err
        = MPI_Allreduce(sendbuf, recvbuf, count, MPI_DOUBLE, op, comm_images_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR(
            "MPI_Allreduce(double*, double*) of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::reduce(double* buf, int count, MPI_Op op, const int root) const
{
    double* recvbuf = new double[count];
    int ret         = reduce(buf, recvbuf, count, op, root);

    if (mype_spin_ == root) memcpy(buf, recvbuf, count * sizeof(double));
    delete[] recvbuf;

    return ret;
}

int MGmol_MPI::reduce(float* buf, int count, MPI_Op op, const int root) const
{
    float* recvbuf = new float[count];
    int ret        = reduce(buf, recvbuf, count, op, root);

    if (mype_spin_ == root) memcpy(buf, recvbuf, count * sizeof(float));
    delete[] recvbuf;

    return ret;
}

int MGmol_MPI::reduce(int* buf, int count, MPI_Op op, const int root) const
{
    int* recvbuf = new int[count];
    int ret      = reduce(buf, recvbuf, count, op, root);

    if (mype_spin_ == root) memcpy(buf, recvbuf, count * sizeof(int));
    delete[] recvbuf;

    return ret;
}

int MGmol_MPI::reduce(short* buf, int count, MPI_Op op, const int root) const
{
    short* recvbuf = new short[count];
    int ret        = reduce(buf, recvbuf, count, op, root);

    if (mype_spin_ == root) memcpy(buf, recvbuf, count * sizeof(short));
    delete[] recvbuf;

    return ret;
}

int MGmol_MPI::allreduce(double* buf, int count, MPI_Op op) const
{
    double* recvbuf = new double[count];
    int ret         = allreduce(buf, recvbuf, count, op);

    memcpy(buf, recvbuf, count * sizeof(double));
    delete[] recvbuf;

    return ret;
}

int MGmol_MPI::allreduce(float* buf, int count, MPI_Op op) const
{
    float* recvbuf = new float[count];
    int ret        = allreduce(buf, recvbuf, count, op);

    memcpy(buf, recvbuf, count * sizeof(float));
    delete[] recvbuf;

    return ret;
}

int MGmol_MPI::allreduce(int* buf, int count, MPI_Op op) const
{
    int* recvbuf = new int[count];
    int ret      = allreduce(buf, recvbuf, count, op);

    memcpy(buf, recvbuf, count * sizeof(int));
    delete[] recvbuf;

    return ret;
}

int MGmol_MPI::allreduceGlobal(int* buf, int count, MPI_Op op) const
{
    int* recvbuf = new int[count];
    int ret      = allreduceGlobal(buf, recvbuf, count, op);

    memcpy(buf, recvbuf, count * sizeof(int));
    delete[] recvbuf;

    return ret;
}

int MGmol_MPI::allreduceGlobal(double* buf, int count, MPI_Op op) const
{
    double* recvbuf = new double[count];
    int ret         = allreduceGlobal(buf, recvbuf, count, op);

    memcpy(buf, recvbuf, count * sizeof(double));
    delete[] recvbuf;

    return ret;
}

int MGmol_MPI::allreduce(short* buf, int count, MPI_Op op) const
{
    short* recvbuf = new short[count];
    int ret        = allreduce(buf, recvbuf, count, op);

    memcpy(buf, recvbuf, count * sizeof(short));
    delete[] recvbuf;

    return ret;
}

int MGmol_MPI::allreduceImages(double* buf, int count, MPI_Op op) const
{
    double* recvbuf = new double[count];
    int ret         = allreduceImages(buf, recvbuf, count, op);

    memcpy(buf, recvbuf, count * sizeof(double));
    delete[] recvbuf;

    return ret;
}

int MGmol_MPI::allreduceSpin(
    double* sendbuf, double* recvbuf, int count, MPI_Op op) const
{
    if (nspin_ == 1) return 0;

    assert(comm_different_spin_ != MPI_COMM_NULL);
    int mpi_err = MPI_Allreduce(
        sendbuf, recvbuf, count, MPI_DOUBLE, op, comm_different_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR(
            "MPI_Allreduce(double*, double*) of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::allreduceSpin(
    float* sendbuf, float* recvbuf, int count, MPI_Op op) const
{
    if (nspin_ == 1) return 0;

    assert(comm_different_spin_ != MPI_COMM_NULL);
    int mpi_err = MPI_Allreduce(
        sendbuf, recvbuf, count, MPI_FLOAT, op, comm_different_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR(
            "MPI_Allreduce(double*, double*) of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::exchangeDataSpin(
    double* localdata, double* remotedata, int count) const
{
    if (nspin_ == 1) return 0;

    assert(localdata != nullptr);
    assert(remotedata != nullptr);
    assert(myspin_ == 0 || myspin_ == 1);

    assert(comm_different_spin_ != MPI_COMM_NULL);
    MPI_Request requestr;
    int src = (myspin_ + 1) % 2;
    int dst = src;
    assert(dst == 0 || dst == 1);
    int mpi_err = MPI_Irecv(remotedata, count, MPI_DOUBLE, src, myspin_,
        comm_different_spin_, &requestr);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MGmol_MPI::exchangeDataSpin(), MPI_Irecv failed!!!");
    }
    MPI_Request requests;
    mpi_err = MPI_Isend(localdata, count, MPI_DOUBLE, dst, (myspin_ + 1) % 2,
        comm_different_spin_, &requests);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MGmol_MPI::exchangeDataSpin(), MPI_Isend failed!!!");
    }

    mpi_err = MPI_Wait(&requestr, MPI_STATUS_IGNORE);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MGmol_MPI::exchangeDataSpin() !!!");
    }
    // barrier();

    return mpi_err;
}

int MGmol_MPI::exchangeDataSpin(
    float* localdata, float* remotedata, int count) const
{
    if (nspin_ == 1) return 0;

    assert(localdata != nullptr);
    assert(remotedata != nullptr);
    assert(myspin_ == 0 || myspin_ == 1);

    assert(comm_different_spin_ != MPI_COMM_NULL);
    MPI_Request requestr;
    int src = (myspin_ + 1) % 2;
    int dst = src;
    assert(dst == 0 || dst == 1);
    int mpi_err = MPI_Irecv(remotedata, count, MPI_FLOAT, src, myspin_,
        comm_different_spin_, &requestr);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("exchangeDataSpin(), MPI_Irecv failed!!!");
    }
    MPI_Request requests;
    mpi_err = MPI_Isend(localdata, count, MPI_FLOAT, dst, (myspin_ + 1) % 2,
        comm_different_spin_, &requests);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("exchangeDataSpin(), MPI_Isend failed!!!");
    }

    mpi_err = MPI_Wait(&requestr, MPI_STATUS_IGNORE);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MGmol_MPI::exchangeDataSpin() !!!");
    }
    // barrier();

    return mpi_err;
}

bool MGmol_MPI::compareSpin(const double val)
{
    double otherval;
    double myval = val;
    exchangeDataSpin(&myval, &otherval, 1);
    if (fabs(val - otherval) > 1.e-8)
    {
        if (instancePE0())
        {
            MGMOL_MPI_ERROR("val=" << val << ", other val=" << otherval);
        }
    }
    return true;
}

int MGmol_MPI::send(double* sendbuf, int count, int dest) const
{
    int mpi_err = MPI_Send(sendbuf, count, MPI_DOUBLE, dest, 0, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Send() of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::recv(double* recvbuf, int count, int src) const
{
    MPI_Status status;
    int mpi_err
        = MPI_Recv(recvbuf, count, MPI_DOUBLE, src, 0, comm_spin_, &status);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Recv() of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::send(float* sendbuf, int count, int dest) const
{
    int mpi_err = MPI_Send(sendbuf, count, MPI_FLOAT, dest, 0, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Send() of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::recv(float* recvbuf, int count, int src) const
{
    MPI_Status status;
    int mpi_err
        = MPI_Recv(recvbuf, count, MPI_FLOAT, src, 0, comm_spin_, &status);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Recv() of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::send(int* sendbuf, int count, int dest) const
{
    int mpi_err = MPI_Send(sendbuf, count, MPI_INT, dest, 0, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Send() of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::recv(int* recvbuf, int count, int src) const
{
    MPI_Status status;
    int mpi_err
        = MPI_Recv(recvbuf, count, MPI_INT, src, 0, comm_spin_, &status);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Recv() of size " << count << "!!!");
    }
    return mpi_err;
}

int MGmol_MPI::bcast(std::string& common_string) const
{
    return bcast(common_string, comm_spin_);
}

int MGmol_MPI::bcastGlobal(std::string& common_string) const
{
    return bcast(common_string, comm_global_);
}

int MGmol_MPI::bcast(std::string& common_string, MPI_Comm comm) const
{
    short size_str = (short)common_string.size();
    int mpi_err    = MPI_Bcast(&size_str, 1, MPI_SHORT, 0, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MGmol_MPI::bcast(), ERROR in bcast of 1 short!!!");
    }
    assert(size_str < 256);

    char* buffer = new char[size_str + 1];
    if (mype_spin_ == 0)
    {
        common_string.copy(buffer, std::string::npos);
        buffer[common_string.length()] = 0;
    }
    mpi_err = MPI_Bcast(buffer, size_str + 1, MPI_CHAR, 0, comm);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Bcast(char*, int) of size " << size_str << "!!!");
    }

    common_string.assign(&buffer[0], size_str);
    delete[] buffer;
    return mpi_err;
}

int MGmol_MPI::gather(
    double* sendbuf, int s_count, double* recvbuf, int recvbufsize, int root)
{
    if (s_count * size_ != recvbufsize)
        MGMOL_MPI_ERROR("incorrect array sizes!");

    int count   = s_count;
    int mpi_err = MPI_Gather(sendbuf, count, MPI_DOUBLE, recvbuf, count,
        MPI_DOUBLE, root, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_gather for int in MGmol_MPI::gather() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::gather(
    float* sendbuf, int s_count, float* recvbuf, int recvbufsize, int root)
{
    if (s_count * size_ != recvbufsize)
        MGMOL_MPI_ERROR("incorrect array sizes!");

    int count   = s_count;
    int mpi_err = MPI_Gather(
        sendbuf, count, MPI_FLOAT, recvbuf, count, MPI_FLOAT, root, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_gather for int in MGmol_MPI::gather() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::gather(
    int* sendbuf, int s_count, int* recvbuf, int recvbufsize, int root)
{
    if (s_count * size_ != recvbufsize)
        MGMOL_MPI_ERROR("incorrect array sizes!");

    int count   = s_count;
    int mpi_err = MPI_Gather(
        sendbuf, count, MPI_INT, recvbuf, count, MPI_INT, root, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_gather for int in MGmol_MPI::gather() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::allGather(
    std::vector<double>& sendbuf, std::vector<double>& recvbuf)
{
    if (sendbuf.size() * size_ != recvbuf.size())
        MGMOL_MPI_ERROR("incorrect array sizes!");
    int count   = (int)sendbuf.size();
    int mpi_err = MPI_Allgather(&sendbuf[0], count, MPI_DOUBLE, &recvbuf[0],
        count, MPI_DOUBLE, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR(
            "MPI_Allgather for double in MGmol_MPI::allGather() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::allGather(std::vector<int>& sendbuf, std::vector<int>& recvbuf)
{
    if (sendbuf.size() * size_ != recvbuf.size())
        MGMOL_MPI_ERROR("incorrect array sizes!");

    int count   = (int)sendbuf.size();
    int mpi_err = MPI_Allgather(
        &sendbuf[0], count, MPI_INT, &recvbuf[0], count, MPI_INT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather for int in MGmol_MPI::allGather() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::allGather(
    double* sendbuf, int s_count, double* recvbuf, int recvbufsize)
{
    if (s_count * size_ != recvbufsize)
        MGMOL_MPI_ERROR("incorrect array sizes!");

    int count   = s_count;
    int mpi_err = MPI_Allgather(
        sendbuf, count, MPI_DOUBLE, recvbuf, count, MPI_DOUBLE, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather for int in MGmol_MPI::allGather() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::allGather(
    float* sendbuf, int s_count, float* recvbuf, int recvbufsize)
{
    if (s_count * size_ != recvbufsize)
        MGMOL_MPI_ERROR("incorrect array sizes!");

    int count   = s_count;
    int mpi_err = MPI_Allgather(
        sendbuf, count, MPI_FLOAT, recvbuf, count, MPI_FLOAT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather for int in MGmol_MPI::allGather() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::allGather(
    int* sendbuf, int s_count, int* recvbuf, int recvbufsize)
{
    if (s_count * size_ != recvbufsize)
        MGMOL_MPI_ERROR("incorrect array sizes!");

    int count   = s_count;
    int mpi_err = MPI_Allgather(
        sendbuf, count, MPI_INT, recvbuf, count, MPI_INT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather for int in MGmol_MPI::allGather() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::allGatherImages(
    std::vector<double>& sendbuf, std::vector<double>& recvbuf)
{
    assert(comm_images_ != MPI_COMM_NULL);
    assert(sendbuf.size() * nimages_ == recvbuf.size());
    int count   = (int)sendbuf.size();
    int mpi_err = MPI_Allgather(&sendbuf[0], count, MPI_DOUBLE, &recvbuf[0],
        count, MPI_DOUBLE, comm_images_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR(
            "MPI_Allgather for double in MGmol_MPI::allGather() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::allGatherImages(
    std::vector<int>& sendbuf, std::vector<int>& recvbuf)
{
    assert(comm_images_ != MPI_COMM_NULL);
    if (sendbuf.size() * nimages_ != recvbuf.size())
    {
        std::cerr << "sendbuf.size()=" << sendbuf.size()
                  << ", recvbuf.size()=" << recvbuf.size()
                  << ", nimages_=" << nimages_ << std::endl;
    }
    assert(sendbuf.size() * nimages_ == recvbuf.size());
    int count   = (int)sendbuf.size();
    int mpi_err = MPI_Allgather(
        &sendbuf[0], count, MPI_INT, &recvbuf[0], count, MPI_INT, comm_images_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather for int in MGmol_MPI::allGather() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::allGatherImages(
    std::vector<short>& sendbuf, std::vector<short>& recvbuf)
{
    assert(comm_images_ != MPI_COMM_NULL);
    if (sendbuf.size() * nimages_ != recvbuf.size())
    {
        std::cerr << "sendbuf.size()=" << sendbuf.size()
                  << ", recvbuf.size()=" << recvbuf.size()
                  << ", nimages_=" << nimages_ << std::endl;
    }
    assert(sendbuf.size() * nimages_ == recvbuf.size());
    int count   = (int)sendbuf.size();
    int mpi_err = MPI_Allgather(&sendbuf[0], count, MPI_SHORT, &recvbuf[0],
        count, MPI_SHORT, comm_images_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather for int in MGmol_MPI::allGather() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::allGatherV(
    std::vector<double>& sendbuf, std::vector<double>& recvbuf)
{
    int sendcount   = (int)sendbuf.size();
    int* recvcounts = new int[size_];
    int mpi_err     = MPI_Allgather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather in MGmol_MPI::allGatherV() !!!");
    }

    int* displs = new int[size_];
    displs[0]   = 0;
    for (int i = 1; i < size_; i++)
    {
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    recvbuf.resize(displs[size_ - 1] + recvcounts[size_ - 1]);
    // assert( recvbuf.size()==displs[size_-1]+recvcounts[size_-1] );

    mpi_err = MPI_Allgatherv(&sendbuf[0], sendcount, MPI_DOUBLE, &recvbuf[0],
        recvcounts, displs, MPI_DOUBLE, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgatherv in MGmol_MPI::allGatherV() !!!");
    }

    delete[] displs;
    delete[] recvcounts;

    return mpi_err;
}

int MGmol_MPI::allGatherV(
    std::vector<float>& sendbuf, std::vector<float>& recvbuf)
{
    int sendcount   = (int)sendbuf.size();
    int* recvcounts = new int[size_];
    int mpi_err     = MPI_Allgather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather in MGmol_MPI::allGatherV() !!!");
    }

    int* displs = new int[size_];
    displs[0]   = 0;
    for (int i = 1; i < size_; i++)
    {
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    recvbuf.resize(displs[size_ - 1] + recvcounts[size_ - 1]);
    // assert( recvbuf.size()==displs[size_-1]+recvcounts[size_-1] );

    mpi_err = MPI_Allgatherv(&sendbuf[0], sendcount, MPI_FLOAT, &recvbuf[0],
        recvcounts, displs, MPI_FLOAT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgatherv in MGmol_MPI::allGatherV() !!!");
    }

    delete[] displs;
    delete[] recvcounts;

    return mpi_err;
}

int MGmol_MPI::allGatherV(std::vector<int>& sendbuf, std::vector<int>& recvbuf)
{
    int sendcount   = (int)sendbuf.size();
    int* recvcounts = new int[size_];
    int mpi_err     = MPI_Allgather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather in MGmol_MPI::allGatherV() !!!");
    }

    int* displs = new int[size_];
    displs[0]   = 0;
    for (int i = 1; i < size_; i++)
    {
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    recvbuf.resize(displs[size_ - 1] + recvcounts[size_ - 1]);
    // assert( recvbuf.size()==displs[size_-1]+recvcounts[size_-1] );

    mpi_err = MPI_Allgatherv(&sendbuf[0], sendcount, MPI_INT, &recvbuf[0],
        recvcounts, displs, MPI_INT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgatherv in MGmol_MPI::allGatherV() !!!");
    }

    delete[] displs;
    delete[] recvcounts;

    return mpi_err;
}

int MGmol_MPI::allGatherV(
    std::vector<short>& sendbuf, std::vector<short>& recvbuf)
{
    int sendcount   = (int)sendbuf.size();
    int* recvcounts = new int[size_];
    int mpi_err     = MPI_Allgather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather in MGmol_MPI::allGatherV() !!!");
    }

    int* displs = new int[size_];
    displs[0]   = 0;
    for (int i = 1; i < size_; i++)
    {
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    recvbuf.resize(displs[size_ - 1] + recvcounts[size_ - 1]);
    // assert( recvbuf.size()==displs[size_-1]+recvcounts[size_-1] );

    mpi_err = MPI_Allgatherv(&sendbuf[0], sendcount, MPI_SHORT, &recvbuf[0],
        recvcounts, displs, MPI_SHORT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgatherv in MGmol_MPI::allGatherV() !!!");
    }

    delete[] displs;
    delete[] recvcounts;

    return mpi_err;
}

int MGmol_MPI::allGatherV(
    std::vector<unsigned short>& sendbuf, std::vector<unsigned short>& recvbuf)
{
    int sendcount   = (int)sendbuf.size();
    int* recvcounts = new int[size_];
    int mpi_err     = MPI_Allgather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather in MGmol_MPI::allGatherV() !!!");
    }

    int* displs = new int[size_];
    displs[0]   = 0;
    for (int i = 1; i < size_; i++)
    {
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    recvbuf.resize(displs[size_ - 1] + recvcounts[size_ - 1]);
    // assert( recvbuf.size()==displs[size_-1]+recvcounts[size_-1] );

    mpi_err = MPI_Allgatherv(&sendbuf[0], sendcount, MPI_UNSIGNED_SHORT,
        &recvbuf[0], recvcounts, displs, MPI_UNSIGNED_SHORT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgatherv in MGmol_MPI::allGatherV() !!!");
    }

    delete[] displs;
    delete[] recvcounts;

    return mpi_err;
}

int MGmol_MPI::allGatherV(
    std::vector<std::string>& sendbuf, std::vector<std::string>& recvbuf)
{
    ////
    int vcount = (int)sendbuf.size();
    allreduce(&vcount, 1, MPI_SUM);

    int sendcount;
    int totchars;
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
    allGatherV(locstrlen, strLen);

    // convert string to char array
    //    vector<char>charStr(sendcount);
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

    int* recvcounts = new int[size_];
    int mpi_err     = MPI_Allgather(
        &sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgather in MGmol_MPI::allGatherV() !!!");
    }

    int* displs = new int[size_];
    displs[0]   = 0;
    totchars    = recvcounts[size_ - 1];
    for (int i = 1; i < size_; i++)
    {
        displs[i] = displs[i - 1] + recvcounts[i - 1];
        totchars += recvcounts[i - 1];
    }

    char* recvdata = new char[totchars];
    mpi_err = MPI_Allgatherv(&charStr[0], sendcount, MPI_CHAR, &recvdata[0],
        recvcounts, displs, MPI_CHAR, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_Allgatherv in MGmol_MPI::allGatherV() !!!");
    }

    // reset recvbuf
    recvbuf.clear();
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

    delete[] displs;
    delete[] recvcounts;
    delete[] charStr;
    delete[] recvdata;

    return mpi_err;
}

int MGmol_MPI::gatherV(
    std::vector<int>& sendbuf, std::vector<int>& recvbuf, const int root)
{
    return mgmol_tools::gatherV(sendbuf, recvbuf, root, comm_spin_);
}

int MGmol_MPI::gatherV(std::vector<std::string>& sendbuf,
    std::vector<std::string>& recvbuf, const int root)
{
    return mgmol_tools::gatherV(sendbuf, recvbuf, root, comm_spin_);
}

int MGmol_MPI::gatherV(
    std::vector<double>& sendbuf, std::vector<double>& recvbuf, const int root)
{
    return mgmol_tools::gatherV(sendbuf, recvbuf, root, comm_spin_);
}

int MGmol_MPI::scatter(
    double* sendbuf, int sendbufsize, double* recvbuf, int count, int root)
{
    if (count * size_ != sendbufsize)
        MGMOL_MPI_ERROR("incorrect scatter size!!!");

    int mpi_err = MPI_Scatter(sendbuf, count, MPI_DOUBLE, recvbuf, count,
        MPI_DOUBLE, root, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_gather for int in MGmol_MPI::scatter() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::scatter(
    float* sendbuf, int sendbufsize, float* recvbuf, int count, int root)
{
    if (count * size_ != sendbufsize)
        MGMOL_MPI_ERROR("incorrect scatter size!!!");

    int mpi_err = MPI_Scatter(
        sendbuf, count, MPI_FLOAT, recvbuf, count, MPI_FLOAT, root, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_gather for int in MGmol_MPI::scatter() !!!");
    }

    return mpi_err;
}

int MGmol_MPI::scatter(
    int* sendbuf, int sendbufsize, int* recvbuf, int count, int root)
{
    if (count * size_ != sendbufsize)
        MGMOL_MPI_ERROR("incorrect scatter size!!!");

    int mpi_err = MPI_Scatter(
        sendbuf, count, MPI_INT, recvbuf, count, MPI_INT, root, comm_spin_);
    if (mpi_err != MPI_SUCCESS)
    {
        MGMOL_MPI_ERROR("MPI_scatter for int in MGmol_MPI::scatter() !!!");
    }

    return mpi_err;
}

void MGmol_MPI::split_allreduce_sums_float(float* array, const int nelements)
{
    split_allreduce_sums_float_tm_.start();

    const int work_size = MIN(MAX_SIZE, nelements);
    float* work_float   = new float[work_size];

    int newsize   = work_size;
    int nnblocks  = nelements / newsize;
    int remainder = (nelements % newsize);

    float* array_ptr = array;

    for (int block = 0; block < nnblocks; block++)
    {

        int rc = MPI_Allreduce(
            array_ptr, work_float, newsize, MPI_FLOAT, MPI_SUM, commSpin());
        if (rc != MPI_SUCCESS)
        {
            MGMOL_MPI_ERROR("MPI_allreduce float sum failed!!!");
        }
        memcpy(array_ptr, work_float, newsize * sizeof(float));

        array_ptr += newsize;
    }

    if (remainder)
    {

        int rc = MPI_Allreduce(
            array_ptr, work_float, remainder, MPI_FLOAT, MPI_SUM, commSpin());
        if (rc != MPI_SUCCESS)
        {
            MGMOL_MPI_ERROR("MPI_allreduce float sum failed!!!");
        }
        memcpy(array_ptr, work_float, remainder * sizeof(float));
    }

    delete[] work_float;

    split_allreduce_sums_float_tm_.stop();
}

void MGmol_MPI::split_allreduce_sums_double(double* array, const int nelements)
{
    split_allreduce_sums_double_tm_.start();

    const int work_size = MIN(MAX_SIZE, nelements);
    double* work_double = new double[work_size];

    int newsize   = work_size;
    int nblocks   = nelements / newsize;
    int remainder = (nelements % newsize);

    double* array_ptr = array;

    for (int block = 0; block < nblocks; block++)
    {

        int rc = MPI_Allreduce(
            array_ptr, work_double, newsize, MPI_DOUBLE, MPI_SUM, commSpin());
        if (rc != MPI_SUCCESS)
        {
            MGMOL_MPI_ERROR("MPI_allreduce double sum failed!!!");
        }
        memcpy(array_ptr, work_double, newsize * sizeof(double));

        array_ptr += newsize;
    }

    if (remainder)
    {

        int rc = MPI_Allreduce(
            array_ptr, work_double, remainder, MPI_DOUBLE, MPI_SUM, commSpin());
        if (rc != MPI_SUCCESS)
        {
            MGMOL_MPI_ERROR("MPI_allreduce double sum failed!!!");
        }
        memcpy(array_ptr, work_double, remainder * sizeof(double));
    }

    delete[] work_double;

    split_allreduce_sums_double_tm_.stop();
}

void MGmol_MPI::split_allreduce_sums_int(int* array, const int nelements)
{
    int* work_int = new int[MAX_SIZE];

    int newsize   = MAX_SIZE;
    int nblocks   = nelements / newsize;
    int remainder = (nelements % newsize);

    int* array_ptr = array;

    for (int block = 0; block < nblocks; block++)
    {

        int rc = MPI_Allreduce(
            array_ptr, work_int, newsize, MPI_INT, MPI_SUM, commSpin());
        if (rc != MPI_SUCCESS)
        {
            MGMOL_MPI_ERROR("MPI_allreduce int sum failed!!!");
        }
        memcpy(array_ptr, work_int, newsize * sizeof(int));

        array_ptr += newsize;
    }

    if (remainder)
    {

        int rc = MPI_Allreduce(
            array_ptr, work_int, remainder, MPI_INT, MPI_SUM, commSpin());
        if (rc != MPI_SUCCESS)
        {
            MGMOL_MPI_ERROR("MPI_allreduce int sum failed!!!");
        }
        memcpy(array_ptr, work_int, remainder * sizeof(int));
    }

    delete[] work_int;
}

void MGmol_MPI::split_allreduce_sums_short(
    short int* array, const int nelements)
{
    assert(nelements > 0);

    short* work_short = new short[MAX_SIZE];

    int newsize   = MAX_SIZE;
    int nblocks   = nelements / newsize;
    int remainder = (nelements % newsize);

    short int* array_ptr = array;

    for (int block = 0; block < nblocks; block++)
    {

        int rc = MPI_Allreduce(
            array_ptr, work_short, newsize, MPI_SHORT, MPI_SUM, commSpin());
        if (rc != MPI_SUCCESS)
        {
            MGMOL_MPI_ERROR("ERROR MGmol_MPI::split_allreduce_sums_short: "
                            "MPI_allreduce int sum failed!!!");
        }
        memcpy(array_ptr, work_short, newsize * sizeof(short));

        array_ptr += newsize;
    }

    if (remainder)
    {

        int rc = MPI_Allreduce(
            array_ptr, work_short, remainder, MPI_SHORT, MPI_SUM, commSpin());
        if (rc != MPI_SUCCESS)
        {
            MGMOL_MPI_ERROR("ERROR MGmol_MPI::split_allreduce_sums_short: "
                            "MPI_allreduce int sum failed!!!");
        }
        memcpy(array_ptr, work_short, remainder * sizeof(short));
    }

    delete[] work_short;
}
