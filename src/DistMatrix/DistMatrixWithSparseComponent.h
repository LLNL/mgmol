// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef DISTMATRIXWithSparseComponent_H
#define DISTMATRIXWithSparseComponent_H

#include "DistMatrix.h"
#include "SparseDistMatrix.h"

namespace dist_matrix{

template <class T>
class DistMatrixWithSparseComponent:public DistMatrix<T>
{
private:
    MPI_Comm comm_;
    
    dist_matrix::SparseDistMatrix<T>* sparse_;
    
    RemoteTasksDistMatrix<T>* rtasks_distmatrix_;
    
    int target_nb_tasks_per_partition_;

    DistMatrixWithSparseComponent& operator=(const DistMatrixWithSparseComponent& mat){}
    DistMatrixWithSparseComponent(const DistMatrixWithSparseComponent& mat){}
public:
    DistMatrixWithSparseComponent(const std::string& name,
                                  const BlacsContext& bc,
                                  const int m, const int n,
                                  MPI_Comm comm,
                                  RemoteTasksDistMatrix<T>* rtasks_distmatrix=NULL,
                                  const int target_nb_tasks_per_partition=0):
        DistMatrix<T>(name,bc,m,n),
        comm_(comm),
        rtasks_distmatrix_(rtasks_distmatrix),
        target_nb_tasks_per_partition_(target_nb_tasks_per_partition)
    {
        sparse_=new SparseDistMatrix<T>(
                   comm, 
                   *this,
                   rtasks_distmatrix,
                   target_nb_tasks_per_partition);
    }
    
    DistMatrixWithSparseComponent(const std::string& name,
                                  const int m, const int n,
                                  MPI_Comm comm,
                                  RemoteTasksDistMatrix<T>* rtasks_distmatrix=NULL,
                                  const int target_nb_tasks_per_partition=0):
        DistMatrix<T>(name,m,n),
        comm_(comm),
        rtasks_distmatrix_(rtasks_distmatrix),
        target_nb_tasks_per_partition_(target_nb_tasks_per_partition)
    {
        sparse_=new SparseDistMatrix<T>(
                   comm, 
                   *this,
                   rtasks_distmatrix,
                   target_nb_tasks_per_partition);
    }

    DistMatrixWithSparseComponent(const std::string& name,
                                  const int m,
                                  MPI_Comm comm):
        DistMatrix<T>(name,m,m),
        comm_(comm)
    {
        sparse_=new SparseDistMatrix<T>(comm, *this);
        
        rtasks_distmatrix_ = sparse_->remoteTasksDistMatrix();
        target_nb_tasks_per_partition_ = SparseDistMatrix<T>::numTasksPerPartitioning();
    }
    
    ~DistMatrixWithSparseComponent()
    {
        delete sparse_;
    }
    
    SparseDistMatrix<T>& sparse(){ return *sparse_; }
    
    DistMatrixWithSparseComponent& operator=(const DistMatrix<T>& mat)
    {
        DistMatrix<T>::operator=(mat);
        return *this;
    }
};

}

#endif
