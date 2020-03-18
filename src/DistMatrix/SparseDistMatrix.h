// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// Storage class for matrix distributed over a set of processors
// Each processor may own some contributions to the matrix denoted
// by a set of indices and values. The global matrix is the sum of
// all the contributions. Indices can be duplicated to represent
// multiple contributions on a single processor.

#ifndef MGMOL_SPARSEDISTMATRIX_H
#define MGMOL_SPARSEDISTMATRIX_H

#include "DistMatrix.h"
#include "MPI_DistMatrixCommunicator.h"
#include "RemoteTasksDistMatrix.h"
#include "SquareSubMatrix.h"
#include "Timer.h"

#include <iostream>
#include <map>
#include <memory>
#include <mpi.h>
#include <vector>

namespace dist_matrix
{

template <class T>
class SparseDistMatrix
{
private:
    static int sparse_distmatrix_nb_partitions_;

    static Timer pSumToDistMatrix_tm_;
    static Timer pSumSendRecv_tm_;
    static Timer assign_tm_;
    static Timer maptoarray_tm_;
    static Timer consolidateArray_tm_;

    static unsigned short consolidation_number_;

    const MPI_Comm comm_global_;
    int mype_;
    int npes_;
    int nprow_; // number of PE row in ScaLapack split
    int ntasks_mat_;

    static std::shared_ptr<RemoteTasksDistMatrix<T>> rtasks_distmatrix_;

    std::shared_ptr<MPI_DistMatrixCommunicator> partition_comm_;
    int ntasks_per_partition_;
    int npartitions_;

    // holds "local" data in sparse format,
    // separated by ScaLapack destination process
    std::vector<std::vector<T>> index_and_val_;
    std::vector<std::map<int, T>> map_val_;

    // holds data for ScaLapack calculations
    DistMatrix<T>& mat_;

    void assign(const int, const T* const);
    void maptoarray();
    void array2map();
    void sendRecvData();
    void parallelSumToDistMatrix1();
    void parallelSumToDistMatrix2();

    int Alltoallv(const void* sendbuf, const int* sendcounts,
        const int* sdispls, void* recvbuf, const int* recvcounts,
        const int* rdispls, MPI_Comm comm);

    int Alltoall(const void* sendbuf, int sendcount, void* recvbuf,
        int recvcount, MPI_Comm comm);

    int Isend(const void* buf, int count, int dest, int tag, MPI_Comm comm,
        MPI_Request* request);

    int Irecv(void* buf, int count, int source, int tag, MPI_Comm comm,
        MPI_Request* request);

public:
    static void setNumTasksPerPartitioning(const int nb)
    {
        sparse_distmatrix_nb_partitions_ = nb;
    }

    static int numTasksPerPartitioning()
    {
        return sparse_distmatrix_nb_partitions_;
    }
    static void printConsolidationNumber(std::ostream& os)
    {
        os << "SparseDistMatrix: consolidation_number_="
           << consolidation_number_ << std::endl;
    }

    static void setConsolidationNumber(unsigned short number)
    {
        consolidation_number_ = number;
    }
    // set rtasks_distmatrix = NULL to create one by default. Also
    // set target_nb_tasks_per_partition to <=0 to use npes by default
    SparseDistMatrix<T>(MPI_Comm comm, DistMatrix<T>& mat,
        const int target_nb_tasks_per_partition);
    SparseDistMatrix<T>(MPI_Comm comm, DistMatrix<T>& mat);
    SparseDistMatrix<T>(const SparseDistMatrix<T>& spdistmat);

    ~SparseDistMatrix<T>();

    void addData(const std::vector<T>& data, const int ld, const int ilow,
        const int ihi, const int jlow, const int jhi);
    void addData(const std::vector<T>& data, const int ld, const int ilow,
        const int ihi, const int jlow, const int jhi,
        const std::vector<int>& gids);
    void addData(const SquareSubMatrix<T>& mat, const double tol = 0.);

    void push_back(const int index1, const int index2, const T val);

    void parallelSumToDistMatrix();

    void print(std::ostream& os) const;

    static void printTimers(std::ostream& os)
    {
        pSumToDistMatrix_tm_.print(os);
        pSumSendRecv_tm_.print(os);
        maptoarray_tm_.print(os);
        assign_tm_.print(os);
        consolidateArray_tm_.print(os);
    }
    void scal(const T alpha);
    void setPartitioning(const int target_nb_tasks_per_partition);
    void consolidateArray();
    size_t size() const;

    void printNpartitions(std::ostream& os)
    {
        os << "SparseDistMatrix: npartitions_=" << npartitions_ << std::endl;
    }

    void clearIndexAndVal()
    {
        int n = (int)index_and_val_.size();
        for (int i = 0; i < n; i++)
        {
            index_and_val_[i].clear();
        }
    }

    void clearData()
    {
        int n = (int)index_and_val_.size();
        for (int i = 0; i < n; i++)
        {
            index_and_val_[i].clear();
        }
        n = (int)map_val_.size();
        for (int i = 0; i < n; i++)
        {
            map_val_[i].clear();
        }
    }
};

template <class T>
int SparseDistMatrix<T>::sparse_distmatrix_nb_partitions_ = 128;

template <class T>
unsigned short SparseDistMatrix<T>::consolidation_number_ = 8;

template <class T>
Timer SparseDistMatrix<T>::pSumToDistMatrix_tm_(
    "SparseDistMatrix::pSumToDistMatrix");
template <class T>
Timer SparseDistMatrix<T>::pSumSendRecv_tm_("SparseDistMatrix::pSumSendRecv");
template <class T>
Timer SparseDistMatrix<T>::assign_tm_("SparseDistMatrix::assign");
template <class T>
Timer SparseDistMatrix<T>::maptoarray_tm_("SparseDistMatrix::maptoarray");
template <class T>
Timer SparseDistMatrix<T>::consolidateArray_tm_(
    "SparseDistMatrix::consolidateArray");

} // namespace

#endif
