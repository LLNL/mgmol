// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_MPI
#define MGMOL_MPI

#include "Timer.h"

#include <mpi.h>

#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>

#ifndef MGMOL_MPI_ERROR

#define MGMOL_MPI_ERROR(X)                                                     \
    do                                                                         \
    {                                                                          \
        std::cerr << "ERROR in file " << __FILE__ << " at line " << __LINE__   \
                  << std::endl;                                                \
        std::cerr << "Error Message: " << X << std::endl;                      \
        abort();                                                               \
    } while (0)

#endif

class MGmol_MPI
{
private:
    static MGmol_MPI* pinstance_;

    // communicator including both spins
    static MPI_Comm comm_global_;

    // communicator for all processes dedicated to a particular spin
    static MPI_Comm comm_spin_;

    // communicator to exchange data between different spin
    static MPI_Comm comm_different_spin_;

    static MPI_Comm comm_images_;

    // global rank
    static int mype_;
    static int mype_spin_;
    static int size_;

    static short myspin_;

    // number of spin (one -- without spin --- or two --- with spin)
    static short nspin_;

    static short myimage_;
    static short nimages_;

    static std::ostream* os_;

    MGmol_MPI();

    static void setupComm(
        const MPI_Comm comm, const bool with_spin, const int nimages);
    static Timer split_allreduce_sums_double_tm_;
    static Timer split_allreduce_sums_float_tm_;

public:
    static MGmol_MPI* instance()
    {
        assert(comm_global_ != MPI_COMM_NULL);
        if (pinstance_ == nullptr)
        {
            pinstance_ = new MGmol_MPI();
        }
        return pinstance_;
    }

    static void deleteInstance()
    {
        delete pinstance_;
        pinstance_ = nullptr;
    }

    static void setup(const MPI_Comm comm, std::ostream& os,
        const bool with_spin = false, const int nimages = 1)
    {
        assert(pinstance_ == 0);

        os_ = &os;

        short wspin = (short)with_spin;
        MPI_Bcast(&wspin, 1, MPI_SHORT, 0, comm);

        setupComm(comm, (bool)wspin, nimages);
    }

    static void printTimers(std::ostream& os);

    MPI_Comm commGlobal() const { return comm_global_; }

    MPI_Comm commSameSpin() const { return comm_spin_; }

    MPI_Comm commSpin() const { return comm_spin_; }

    int mypeSpin() const { return mype_spin_; }
    int mypeGlobal() const { return mype_; }
    int myspin() const { return myspin_; }
    int myimage() const
    {
        assert(myimage_ < nimages_);
        assert(myimage_ >= 0);
        return myimage_;
    }
    int nspin() const { return nspin_; }
    bool instancePE0() const { return (mype_spin_ == 0); }
    bool PE0() const { return (mype_ == 0); }

    int size() const { return size_; }

    void barrier() const { MPI_Barrier(comm_global_); }
    int bcast(std::string& common_string, MPI_Comm comm) const;
    int bcast(short* val, int size = 1, int root = 0) const;
    int bcast(int* val, int size = 1, int root = 0) const;
    int bcast(double* val, int size = 1, int root = 0) const;
    int bcast(float* val, int size = 1, int root = 0) const;
    int bcast(char* val, int size = 1, int root = 0) const;
    int bcast(std::string&) const;
    int bcastGlobal(short* val, int size = 1, int root = 0) const;
    int bcastGlobal(int* val, int size = 1, int root = 0) const;
    int bcastGlobal(double* val, int size = 1, int root = 0) const;
    int bcastGlobal(std::string&) const;

    int getRankMinVal(const double val, int* rank);
    int getRankMaxVal(const float val, int* rank);

    int reduce(short* sendbuf, short* recvbuf, int count, MPI_Op op,
        const int root) const;
    int reduce(
        int* sendbuf, int* recvbuf, int count, MPI_Op op, const int root) const;
    int reduce(double* sendbuf, double* recvbuf, int count, MPI_Op op,
        const int root) const;
    int reduce(float* sendbuf, float* recvbuf, int count, MPI_Op op,
        const int root) const;

    int allreduce(short* sendbuf, short* recvbuf, int count, MPI_Op op) const;
    int allreduce(int* sendbuf, int* recvbuf, int count, MPI_Op op) const;
    int allreduce(double* sendbuf, double* recvbuf, int count, MPI_Op op) const;
    int allreduce(float* sendbuf, float* recvbuf, int count, MPI_Op op) const;
    int allreduceGlobal(int* sendbuf, int* recvbuf, int count, MPI_Op op) const;
    int allreduceGlobal(
        double* sendbuf, double* recvbuf, int count, MPI_Op op) const;

    int allreduce(short* buf, int count, MPI_Op op) const;
    int allreduce(int* buf, int count, MPI_Op op) const;
    int allreduce(double* buf, int count, MPI_Op op) const;
    int allreduce(float* buf, int count, MPI_Op op) const;
    int allreduceGlobal(double* buf, int count, MPI_Op op) const;
    int allreduceGlobal(int* buf, int count, MPI_Op op) const;

    int reduce(short* buf, int count, MPI_Op op, const int root) const;
    int reduce(int* buf, int count, MPI_Op op, const int root) const;
    int reduce(double* buf, int count, MPI_Op op, const int root) const;
    int reduce(float* buf, int count, MPI_Op op, const int root) const;

    int allreduceSpin(
        double* sendbuf, double* recvbuf, int count, MPI_Op op) const;
    int allreduceSpin(
        float* sendbuf, float* recvbuf, int count, MPI_Op op) const;
    int exchangeDataSpin(
        double* localdata, double* remotedata, int count) const;
    int exchangeDataSpin(float* localdata, float* remotedata, int count) const;
    int send(double* sendbuf, int count, int dest) const;
    int recv(double* recvbuf, int count, int src) const;
    int send(float* sendbuf, int count, int dest) const;
    int recv(float* recvbuf, int count, int src) const;
    int send(int* sendbuf, int count, int dest) const;
    int recv(int* recvbuf, int count, int src) const;
    bool compareSpin(const double val);

    int gather(double* sendbuf, int scount, double* recvbuf, int recvbufsize,
        int root = 0);
    int gather(float* sendbuf, int scount, float* recvbuf, int recvbufsize,
        int root = 0);
    int gather(
        int* sendbuf, int scount, int* recvbuf, int recvbufsize, int root = 0);

    int gatherV(
        std::vector<int>& sendbuf, std::vector<int>& recvbuf, const int root);
    int gatherV(std::vector<double>& sendbuf, std::vector<double>& recvbuf,
        const int root);

    int allGatherV(std::vector<double>& sendbuf, std::vector<double>& recvbuf);
    int allGatherV(std::vector<float>& sendbuf, std::vector<float>& recvbuf);
    int allGatherV(std::vector<int>& sendbuf, std::vector<int>& recvbuf);
    int allGatherV(std::vector<short>& sendbuf, std::vector<short>& recvbuf);
    int allGatherV(std::vector<unsigned short>& sendbuf,
        std::vector<unsigned short>& recvbuf);

    int gatherV(std::vector<std::string>& sendbuf,
        std::vector<std::string>& recvbuf, const int root);
    int allGatherV(
        std::vector<std::string>& sendbuf, std::vector<std::string>& recvbuf);

    int allGather(std::vector<double>& sendbuf, std::vector<double>& recvbuf);
    int allGather(std::vector<int>& sendbuf, std::vector<int>& recvbuf);
    int allGather(
        double* sendbuf, int send_count, double* recvbuf, int recvbufsize);
    int allGather(
        float* sendbuf, int send_count, float* recvbuf, int recvbufsize);
    int allGather(int* sendbuf, int send_count, int* recvbuf, int recvbufsize);

    int allGatherImages(
        std::vector<short>& sendbuf, std::vector<short>& recvbuf);
    int allGatherImages(std::vector<int>& sendbuf, std::vector<int>& recvbuf);
    int allGatherImages(
        std::vector<double>& sendbuf, std::vector<double>& recvbuf);

    int scatter(double* sendbuf, int sendbufsize, double* recvbuf, int count,
        int root = 0);
    int scatter(float* sendbuf, int sendbufsize, float* recvbuf, int count,
        int root = 0);
    int scatter(
        int* sendbuf, int sendbufsize, int* recvbuf, int count, int root = 0);

    int allreduceImages(double* buf, int count, MPI_Op op) const;
    int allreduceImages(
        double* sendbuf, double* recvbuf, int count, MPI_Op op) const;

    void split_allreduce_sums_double(double*, const int);
    void split_allreduce_sums_float(float*, const int);
    void split_allreduce_sums_int(int*, const int);
    void split_allreduce_sums_short(short int*, const int);

    void abort() const { MPI_Abort(comm_global_, EXIT_FAILURE); }
};

#endif
