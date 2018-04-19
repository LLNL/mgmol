// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: SparseDistMatrix.h 17 2012-08-15 17:11:26Z jeanluc $
// Storage class for matrix distributed over a set of processors
// Each processor may own some contributions to the matrix denoted
// by a set of indices and values. The global matrix is the sum of
// all the contributions. Indices can be duplicated to represent 
// multiple contributions on a single processor.

#ifndef SPARSEDISTMATRIX_H
#define SPARSEDISTMATRIX_H

#if USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include "DistMatrix.h"
#include "Timer.h"
#include "RemoteTasksDistMatrix.h"
#include "MPI_DistMatrixCommunicator.h"

#include <vector>
#include <map>
#include <iostream>
#include <boost/smart_ptr.hpp>

namespace dist_matrix{

template <class T>
class SparseDistMatrix
{
private:
    static int sparse_distmatrix_nb_partitions_;
    static RemoteTasksDistMatrix<T>** def_rtasks_DistMatrix_ptr_;

    static Timer   pSumToDistMatrix_tm_;
    static Timer   pSumSendRecv_tm_;
    static Timer   assign_tm_;
    static Timer   maptoarray_tm_;
    static Timer   consolidateArray_tm_;
  
    static unsigned short consolidation_number_;

    const MPI_Comm comm_global_;
    int        mype_;
    int        npes_;
    int        nprow_; // number of PE row in ScaLapack split
    int        ntasks_mat_;
    
    bool own_rtasks_distmatrix_;
    RemoteTasksDistMatrix<T>* rtasks_distmatrix_;
    //vector<short> other_tasks_indexes_;
    
    boost::shared_ptr<MPI_DistMatrixCommunicator>   partition_comm_;
    int        ntasks_per_partition_;
    int        npartitions_;
    
    // holds "local" data in sparse format,
    // separated by ScaLapack destination process 
    std::vector<std::vector<T> >  index_and_val_;
    std::vector<std::map<int,T> > map_val_;
    
    // holds data for ScaLapack calculations
    DistMatrix<T>& mat_;

    void assign(const int, const T* const);
    void maptoarray();
    void array2map();
    void sendRecvData();
    void parallelSumToDistMatrix1();
    void parallelSumToDistMatrix2();

public:
    static void setNumTasksPerPartitioning(const int nb)
    {
        sparse_distmatrix_nb_partitions_=nb;
    }
    
    static void setRemoteTasksDistMatrixPtr(RemoteTasksDistMatrix<T>** rtasksDistMat)
    {
       def_rtasks_DistMatrix_ptr_=rtasksDistMat;
    }
    static int numTasksPerPartitioning()
    {
       return sparse_distmatrix_nb_partitions_;
    }
    static RemoteTasksDistMatrix<T>* defaultRemoteTasksDistMatrix()
    {
       return *def_rtasks_DistMatrix_ptr_; 
    }
    static void printConsolidationNumber(std::ostream& os)
    {
       os<<"SparseDistMatrix: consolidation_number_="
         <<consolidation_number_<<endl;
    }
 
    static void setConsolidationNumber(unsigned short number)
    {
       consolidation_number_=number;
    }
    RemoteTasksDistMatrix<T>* remoteTasksDistMatrix()
    {
       return rtasks_distmatrix_; 
    }
    // set rtasks_distmatrix = NULL to create one by default. Also 
    // set target_nb_tasks_per_partition to <=0 to use npes by default
    SparseDistMatrix<T>(MPI_Comm comm, DistMatrix<T>& mat,
                        RemoteTasksDistMatrix<T>* rtasks_distmatrix,
                        const int target_nb_tasks_per_partition);
    SparseDistMatrix<T>(MPI_Comm comm, DistMatrix<T>& mat);
    SparseDistMatrix<T>(const SparseDistMatrix<T>& spdistmat);
    
    ~SparseDistMatrix<T>();

    void push_back(const int index1, const int index2, const T val);
    
    void parallelSumToDistMatrix();
    
    void print(std::ostream& os)const;
    
    static void printTimers(std::ostream& os)
    {
      pSumToDistMatrix_tm_.print(os);
      pSumSendRecv_tm_.print(os);
      maptoarray_tm_.print(os);
      assign_tm_.print(os);
      consolidateArray_tm_.print(os);
    }
    void scal(const double alpha);
    void setPartitioning(const int target_nb_tasks_per_partition);
    void consolidateArray();
    size_t size()const;
    
    void printNpartitions(std::ostream& os)
    {
       os<<"SparseDistMatrix: npartitions_="<<npartitions_<<endl;
    }
 
    void clearIndexAndVal()
    {
      int n=(int)index_and_val_.size();
      for(int i=0;i<n;i++)
      {
        index_and_val_[i].clear();
      }
    }
    
    void clearData()
    {
      int n=(int)index_and_val_.size();
      for(int i=0;i<n;i++)
      {
        index_and_val_[i].clear();
      }
      n=(int)map_val_.size();
      for(int i=0;i<n;i++)
      {
        map_val_[i].clear();
      }
    }
};

} // namespace

#endif
