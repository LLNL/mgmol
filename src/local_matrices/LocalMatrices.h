// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_LOCALMATRICES_H
#define MGMOL_LOCALMATRICES_H

#include "LocalVector.h"
#include "MGmol_blas1.h"

#include "Timer.h"
#include "memory_space.h"
#include "mputils.h"

#include "../global.h"

#ifdef HAVE_BML
extern "C"
{
#include <bml.h>
}
#endif

#include <cassert>
#include <iostream>
#include <vector>

const double tol_mat_elements = 1.e-14;

/* LOCALMATRICES class - matrix entries are accessed in column-major order */

template <typename DataType, typename MemorySpaceType>
class LocalMatrices
{
    int storage_size_;
    std::vector<DataType*> ptr_matrices_;

    void set2zero();
    void allocate();

protected:
    // number of rows in matrices
    const int m_;
    // number of columns in matrices
    const int n_;
    // number of matrices
    const short nmat_;

    std::unique_ptr<DataType, void (*)(DataType*)> storage_;

public:
    LocalMatrices(const short nmat, const int m, const int n);
    LocalMatrices(const LocalMatrices&);

    template <class DataType2>
    void copy(const LocalMatrices<DataType2, MemorySpaceType>& mat);
#ifdef HAVE_BML
    void copy(const bml_matrix_t* A);
#endif

    virtual ~LocalMatrices() {};

    short nmat() const { return nmat_; }

    int n() const { return n_; }

    int m() const { return m_; }

    DataType* getSubMatrix(const int iloc = 0) const
    {
        assert(iloc < (int)ptr_matrices_.size());
        assert(ptr_matrices_[iloc] != NULL);
        return ptr_matrices_[iloc];
    }

    DataType* getRawPtr(const int iloc = 0)
    {
        assert(iloc < (int)ptr_matrices_.size());
        assert(ptr_matrices_[iloc] != NULL);
        return ptr_matrices_[iloc];
    }

    void setValues(DataType* values, const int ld, const int iloc = 0);

#ifdef HAVE_MAGMA
    void assign(const LocalMatrices<DataType, MemorySpace::Host>& src);
    void assign(const LocalMatrices<DataType, MemorySpace::Device>& src);
#endif

    DataType getVal(const int i, const int j, const int iloc = 0)
    {
        assert(i < m_);
        assert(j < n_);
        return ptr_matrices_[iloc][m_ * j + i];
    }

    void scal(const double alpha);
    void axpy(const double alpha, const LocalMatrices& matA)
    {
        Taxpy(storage_size_, alpha, matA.storage_.get(), storage_.get());
    }

    void syrk(const int iloc, const int m, const float* const a, const int lda);
    void syrk(
        const int iloc, const int m, const double* const a, const int lda);
    void gemm(const int iloc, const int ma, const float* const a, const int lda,
        const float* const b, const int ldb);
    void gemm(const int iloc, const int ma, const double* const a,
        const int lda, const double* const b, const int ldb);
    void gemm(const char transa, const char transb, const double alpha,
        const LocalMatrices& matA, const LocalMatrices& matB,
        const double beta);
    void reset()
    {
        memset(storage_.get(), 0, storage_size_ * sizeof(DataType));
    }

    void setValues(const DataType val)
    {
        for (int iloc = 0; iloc < nmat_; iloc++)
        {
            DataType* ssiloc = ptr_matrices_[iloc];
            for (int i = 0; i < m_; i++)
            {
                for (int j = 0; j < n_; j++)
                {
                    ssiloc[i + j * m_] = val;
                }
            }
        }
    }

    void print(std::ostream& os, const int iloc) const
    {
        os << "LocalMatrices for iloc=" << iloc << std::endl;
        os << std::scientific;
        const DataType* const ssiloc = ptr_matrices_[iloc];
        for (int i = 0; i < m_; i++)
        {
            for (int j = 0; j < n_; j++)
            {
                os << ssiloc[i + j * m_] << "\t";
            }
            os << std::endl;
        }
    }

    void printBlock(std::ostream& os, const int iloc, const short bsize) const
    {
        os << "LocalMatrices for iloc=" << iloc << std::endl;
        os << std::scientific;
        const DataType* const ssiloc = ptr_matrices_[iloc];
        for (int i = 0; i < bsize; i++)
        {
            for (int j = 0; j < bsize; j++)
            {
                os << ssiloc[i + j * m_] << "\t";
            }
            os << std::endl;
        }
    }

    void print(std::ostream& os) const
    {
        for (int iloc = 0; iloc < nmat_; iloc++)
            print(os, iloc);
    }

    static void printTimers(std::ostream& /*os*/) { }

    void applyMask(const LocalMatrices& mask);

    void setMaskThreshold(
        const DataType min_threshold, const DataType max_threshold);
    void printBlock(std::ostream& os, const int blocksize);

    void matvec(const LocalVector<DataType, MemorySpaceType>& u,
        LocalVector<DataType, MemorySpaceType>& f, const int iloc = 0);
};

#endif
