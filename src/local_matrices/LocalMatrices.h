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

template <class T>
class LocalMatrices
{
    T* storage_;
    int storage_size_;
    std::vector<T*> ptr_matrices_;

protected:
    const int m_;
    const int n_;
    const short subdiv_;

public:
    LocalMatrices(const short subdiv, const int m, const int n);
    LocalMatrices(const LocalMatrices&);

    template <class T2>
    void copy(const LocalMatrices<T2>& mat);
#ifdef HAVE_BML
    void copy(const bml_matrix_t* A);
#endif

    virtual ~LocalMatrices()
    {
        if (storage_ != nullptr)
        {
            delete[] storage_;
            storage_ = nullptr;
        }
    };

    short subdiv() const { return subdiv_; }

    int n() const { return n_; }

    int m() const { return m_; }

    const T* getSubMatrix(const int iloc = 0) const
    {
        assert(iloc < (int)ptr_matrices_.size());
        assert(ptr_matrices_[iloc] != NULL);
        return ptr_matrices_[iloc];
    }

    T* getSubMatrix(const int iloc = 0)
    {
        assert(iloc < (int)ptr_matrices_.size());
        assert(ptr_matrices_[iloc] != NULL);
        return ptr_matrices_[iloc];
    }

    void setVal(const int i, const int j, const T val, const int iloc = 0)
    {
        assert(i < m_);
        assert(j < n_);
        ptr_matrices_[iloc][m_ * j + i] = val;
    }

    T getVal(const int i, const int j, const int iloc = 0)
    {
        assert(i < m_);
        assert(j < n_);
        return ptr_matrices_[iloc][m_ * j + i];
    }

    void scal(const double alpha) { Tscal(storage_size_, alpha, storage_); }
    void axpy(const double alpha, const LocalMatrices& matA)
    {
        Taxpy(storage_size_, alpha, matA.storage_, storage_);
    }

    template <typename MemorySpaceType>
    void syrk(const int iloc, const int m, const float* const a, const int lda);
    template <typename MemorySpaceType>
    void syrk(
        const int iloc, const int m, const double* const a, const int lda);
    void gemm(const int iloc, const int ma, const float* const a, const int lda,
        const float* const b, const int ldb);
    void gemm(const int iloc, const int ma, const double* const a,
        const int lda, const double* const b, const int ldb);
    void gemm(const char transa, const char transb, const double alpha,
        const LocalMatrices& matA, const LocalMatrices& matB,
        const double beta);
    void reset() { memset(storage_, 0, storage_size_ * sizeof(T)); }

    void setValues(const T val)
    {
        for (int iloc = 0; iloc < subdiv_; iloc++)
        {
            T* ssiloc = ptr_matrices_[iloc];
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
        const T* const ssiloc = ptr_matrices_[iloc];
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
        const T* const ssiloc = ptr_matrices_[iloc];
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
        for (int iloc = 0; iloc < subdiv_; iloc++)
            print(os, iloc);
    }

    static void printTimers(std::ostream& /*os*/) {}

    void applyMask(const LocalMatrices& mask);

    void setMaskThreshold(const T min_threshold, const T max_threshold);
    void printBlock(std::ostream& os, const int blocksize);

    void matvec(const LocalVector<T>& u, LocalVector<T>& f, const int iloc = 0);
};

#endif
