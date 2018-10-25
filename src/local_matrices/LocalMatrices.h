// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_LOCALMATRICES_H
#define MGMOL_LOCALMATRICES_H

#include "MGmol_blas1.h"

#include "SubMatrices.h"
#include "SparseDistMatrix.h"
#include "Timer.h"
#include "mputils.h"

#include "../global.h"

#include <iostream>
#include <vector>
#include <cassert>


const double tol_mat_elements=1.e-14;

/* LOCALMATRICES class - matrix entries are accessed in column-major order */

template <class T>
class LocalMatrices
{
    static Timer fill_dist_matrix_tm_;
    static dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>* remote_tasks_DistMatrix_;
    static short sparse_distmatrix_nb_tasks_per_partitions_;

    T* storage_;
    const int     n_;
    int     storage_size_;
    std::vector<T*> ptr_matrices_;
    
protected:
    const int   m_;
    const short subdiv_;

public:

    LocalMatrices(const short subdiv, const int m, const int n);
    LocalMatrices(const LocalMatrices&);
    
    template <class T2>
    void copy(const LocalMatrices<T2>& mat);

    virtual ~LocalMatrices()
    {
        if( storage_!=0 )
        {
            delete[] storage_;
            storage_=0;
        }
    };
    
    static void registerRemoteTasksDistMatrix(dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>* remote_tasks_DistMatrix)
    {
        assert( remote_tasks_DistMatrix!=0 );
        remote_tasks_DistMatrix_ = remote_tasks_DistMatrix;
    }
    
    short subdiv()const
    {
        return subdiv_;
    }
    
    int n()const
    {
        return n_;
    }
    
    int m()const
    {
        return m_;
    }
    
    const T* getSubMatrix(const int iloc)const
    {
        assert( iloc<(int)ptr_matrices_.size() );
        assert( ptr_matrices_[iloc]!=NULL );
        return ptr_matrices_[iloc];
    }
    
    T* getSubMatrix(const int iloc)
    {
        assert( iloc<(int)ptr_matrices_.size() );
        assert( ptr_matrices_[iloc]!=NULL );
        return ptr_matrices_[iloc];
    }
    
    void addVal(const int iloc, const int index, const T val)
    {
        ptr_matrices_[iloc][index]+=val;
    }

    void addVal(const int iloc, const int i, const int j, 
                const T val)
    {
        assert( i<m_ );
        assert( j<n_ );
        ptr_matrices_[iloc][m_*j+i]+=val;
    }
    
    // use fortran convention to be compatible with BLAS
    T getVal(const int iloc, const int i, const int j)const
    {
        assert( i<m_ );
        assert( j<n_ );
        return ptr_matrices_[iloc][m_*j+i];
    }
    
    void setVal(const int iloc, const int i, const int j, 
                const T val)
    {
        assert( i<m_ );
        assert( j<n_ );
        ptr_matrices_[iloc][m_*j+i]=val;
    }
    
    void setVal2zero(const int iloc, const int i, const int j)
    {
        assert( i<m_ );
        assert( j<n_ );
        ptr_matrices_[iloc][m_*j+i]=0.;
    }
    
    void scal(const double alpha)
    {
        Tscal(storage_size_, alpha, storage_);
    }
    
    void syrk(const int iloc, const int m,
          const float* const a, const int lda);
    void syrk(const int iloc, const int m,
          const double* const a, const int lda);          
    void gemm(const int iloc, const int ma,
          const float* const a, const int lda,
          const float* const b, const int ldb);
    void  gemm(const int iloc, const int ma,
          const double* const a, const int lda,
          const double* const b, const int ldb);          
    void gemm(const char transa, const char transb, 
          const double alpha,const LocalMatrices& matA, 
          const LocalMatrices& matB, const double beta); 
    void reset()
    {
        memset(storage_,0,storage_size_*sizeof(T));
    }
    
    void setValues(const T val)
    {
        for(int iloc=0;iloc<subdiv_;iloc++)
        {
            T* ssiloc=ptr_matrices_[iloc]; 
            for(int i = 0; i < m_; i++)
            { 
                for(int j = 0; j < n_; j++)
                {
                    ssiloc[i + j*m_]=val;
                }
            }
        }
    }
    
    void print(std::ostream& os, const int iloc)const
    {
        os<<"LocalMatrices for iloc="<<iloc<<std::endl;
        os<<scientific;
        const T* const ssiloc=ptr_matrices_[iloc]; 
        for(int i = 0; i < m_; i++)
        { 
            for(int j = 0; j < n_; j++)
            {
                os<<ssiloc[i + j*m_]<<"\t";
            }
            os<<endl;
        }
    }

    void printBlock(std::ostream& os, const int iloc, const short bsize)const
    {
        os<<"LocalMatrices for iloc="<<iloc<<std::endl;
        os<<scientific;
        const T* const ssiloc=ptr_matrices_[iloc]; 
        for(int i = 0; i < bsize; i++)
        { 
            for(int j = 0; j < bsize; j++)
            {
                os<<ssiloc[i + j*m_]<<"\t";
            }
            os<<endl;
        }
    }

    void print(std::ostream& os)const
    {
        for(int iloc=0;iloc<subdiv_;iloc++)
            print(os,iloc);
    }
    
    static void printTimers(std::ostream& os)
    {
        fill_dist_matrix_tm_.print(os);
    }

    void fillSparseDistMatrix(dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sm,
                              const std::vector<std::vector<int> >& global_indexes,
                              const int numst,
                              const double tol=tol_mat_elements)const;
    void fillDistMatrix(dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
                        const std::vector<std::vector<int> >& global_indexes,
                        const double tol=tol_mat_elements)const;
    void init(const dist_matrix::DistMatrix<DISTMATDTYPE>& mat, 
              const std::vector<std::vector<int> >& global_indexes,
              const dist_matrix::SubMatricesIndexing<DISTMATDTYPE>& submat_indexing);
    void applyMask(const LocalMatrices& mask);
    
    void setMaskThreshold(const T min_threshold,
                          const T max_threshold);
};

template <class T>
Timer LocalMatrices<T>::fill_dist_matrix_tm_("LocalMatrices::fill_dist_matrix");

template <class T>
dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>*
    LocalMatrices<T>::remote_tasks_DistMatrix_=0;

template <class T>
short LocalMatrices<T>::sparse_distmatrix_nb_tasks_per_partitions_=256;

#endif

