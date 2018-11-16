// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "LocalMatrices.h"
#include "MGmol_MPI.h"
#include "blas3_c.h"

template <class T>
LocalMatrices<T>::LocalMatrices(const short subdiv, const int m, const int n)
    : subdiv_(subdiv), m_(m), n_(n)
{
    assert(m >= 0);
    assert(n >= 0);
    assert(subdiv >= 1);

    ptr_matrices_.resize(subdiv_);
    storage_size_ = subdiv_ * m_ * n_;

    if (storage_size_ > 0)
    {
        storage_ = new T[storage_size_];
        memset(storage_, 0, storage_size_ * sizeof(T));
        const int dim2 = m_ * n_;
        for (int i = 0; i < subdiv; i++)
            ptr_matrices_[i] = storage_ + dim2 * i;
    }
    else
    {
        storage_ = 0;
    }
}

template <class T>
LocalMatrices<T>::LocalMatrices(const LocalMatrices& mat)
    : subdiv_(mat.subdiv_), m_(mat.m_), n_(mat.n_)
{
    storage_size_ = mat.storage_size_;

    ptr_matrices_.resize(subdiv_);

    if (storage_size_ > 0)
    {
        storage_ = new T[storage_size_];
        memcpy(storage_, mat.storage_, storage_size_ * sizeof(T));
        const int dim2 = m_ * n_;
        for (int i = 0; i < subdiv_; i++)
            ptr_matrices_[i] = storage_ + dim2 * i;
    }
    else
    {
        storage_ = 0;
    }
}

template <class T1>
template <class T2>
void LocalMatrices<T1>::copy(const LocalMatrices<T2>& mat)
{
    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        T1* dst_mat             = ptr_matrices_[iloc];
        const T2* const src_mat = mat.getSubMatrix(iloc);
        const int nelements     = n_ * m_;
        for (int j = 0; j < nelements; j++)
        {
            dst_mat[j] = src_mat[j];
        }
    }
}

// perform the symmetric operation
// C := A'*A
// This one is for debug purposes, and can be removed if unused
template <class T>
void LocalMatrices<T>::syrk(
    const int iloc, const int k, const double* const a, const int lda)
{
    assert(iloc < subdiv_);
    assert(iloc < (int)ptr_matrices_.size());
    assert(k <= lda);
    assert(ptr_matrices_[iloc] != 0);
    assert(m_ > 0);
    assert(m_ == n_);

    T* const ssiloc = ptr_matrices_[iloc];
    assert(ssiloc != 0);

    const double zero = 0.;
    const double one  = 1.;
    const char uplo   = 'l'; // fill lower triangular part
    const char trans  = 't';
    MPsyrk(uplo, trans, m_, k, one, a, lda, zero, ssiloc, m_);
    //    dsyrk(&uplo, &trans, &m_, &k, &one, a, &lda,
    //          &zero, ssiloc, &m_);
}

// perform the symmetric operation
// C := A'*A
template <class T>
void LocalMatrices<T>::syrk(
    const int iloc, const int k, const float* const a, const int lda)
{
    assert(iloc < subdiv_);
    assert(iloc < (int)ptr_matrices_.size());
    assert(k <= lda);
    assert(ptr_matrices_[iloc] != 0);
    assert(m_ > 0);
    assert(m_ == n_);

    T* const ssiloc = ptr_matrices_[iloc];
    assert(ssiloc != 0);

    const double zero = 0.;
    const double one  = 1.;
    const char uplo   = 'l'; // fill lower triangular part
    const char trans  = 't';
    MPsyrk(uplo, trans, m_, k, one, a, lda, zero, ssiloc, m_);
    //    dsyrk(&uplo, &trans, &m_, &k, &one, a, &lda,
    //          &zero, ssiloc, &m_);
}

// This one is for debug purposes, and can be removed if unused
template <class T>
void LocalMatrices<T>::gemm(const int iloc, const int ma, const double* const a,
    const int lda, const double* const b, const int ldb)
{
    assert(iloc < subdiv_);
    assert(iloc < (int)ptr_matrices_.size());
    assert(ma <= lda);
    assert(ptr_matrices_[iloc] != 0);

    T* const c = ptr_matrices_[iloc];
    assert(c != 0);

    MPgemm('t', 'n', m_, n_, ma, 1., a, lda, b, ldb, 0., c, m_);
}

template <class T>
void LocalMatrices<T>::gemm(const int iloc, const int ma, const float* const a,
    const int lda, const float* const b, const int ldb)
{
    assert(iloc < subdiv_);
    assert(iloc < (int)ptr_matrices_.size());
    assert(ma <= lda);
    assert(ptr_matrices_[iloc] != 0);

    T* const c = ptr_matrices_[iloc];
    assert(c != 0);

    MPgemm('t', 'n', m_, n_, ma, 1., a, lda, b, ldb, 0., c, m_);
}

// matrix multiplication
// this = alpha*op(A)*op(B)+beta*this
template <class T>
void LocalMatrices<T>::gemm(const char transa, const char transb,
    const double alpha, const LocalMatrices& matA, const LocalMatrices& matB,
    const double beta)
{
    assert(matA.n() == matB.m());

    const int lda = matA.m();
    const int ldb = matB.m();
    const int nca = matA.n();

    // loop over subdomains
    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        // get matrix data from matA and MatB
        const T* const amat = matA.getSubMatrix(iloc);
        const T* const bmat = matB.getSubMatrix(iloc);
        // get pointer to local storage
        T* const c = ptr_matrices_[iloc];
        assert(c != 0);

        // do matrix multiplication
        MPgemm(transa, transb, m_, n_, nca, alpha, amat, lda, bmat, ldb, beta,
            c, m_);
    }
}

template <class T>
void LocalMatrices<T>::fillSparseDistMatrix(
    dist_matrix::SparseDistMatrix<DISTMATDTYPE>& dst,
    const vector<vector<int>>& global_indexes, const int numst,
    const double tol) const
{
    int* ist = new int[2 * subdiv_];
    T* val   = new T[subdiv_];

    const short chromatic_number = (short)global_indexes[0].size();

    // double loop over colors
    for (short icolor = 0; icolor < chromatic_number; icolor++)
    {
        for (short jcolor = 0; jcolor < chromatic_number; jcolor++)
        {
            int pst_old  = -1;
            int valindex = -1;

            // loop over subdomains
            for (short iloc = 0; iloc < subdiv_; iloc++)
            {
                const int st1 = global_indexes[iloc][icolor];
                //(*MPIdata::sout)<<"icolor="<<icolor<<", pst1="<<pst1<<endl;

                if (st1 != -1)
                {
                    const int st2 = global_indexes[iloc][jcolor];
                    if (st2 != -1)
                    {
                        // unique id for current pair
                        const int pst         = st2 * numst + st1;
                        const T* const ssiloc = getSubMatrix(iloc);
                        const T tmp
                            = ssiloc[jcolor * chromatic_number + icolor];

                        // accumulate values
                        if (fabs(tmp) > tol)
                        {
                            if (pst == pst_old)
                            {
                                assert(valindex >= 0);
                                val[valindex] += tmp;
                            }
                            else
                            {
                                valindex++;
                                val[valindex]         = tmp;
                                ist[2 * valindex]     = st1;
                                ist[2 * valindex + 1] = st2;
                                pst_old               = pst;
                            }
                        }
                    }
                }

            } // iloc

            // assign values
            for (int i = 0; i <= valindex; i++)
                dst.push_back(ist[2 * i], ist[2 * i + 1], (double)val[i]);

        } // jcolor
    } // icolor

    delete[] val;
    delete[] ist;
}

// Add matrix elements corresponding to subdomains at their right place
// into distributed matrix
// Important Note: Neglect contributions smaller than tol!
// (may lead to results dependent on number of CPUs)
template <class T>
void LocalMatrices<T>::fillDistMatrix(
    dist_matrix::DistMatrix<DISTMATDTYPE>& dst,
    const vector<vector<int>>& global_indexes, const double tol) const
{
    // if(onpe0)(*MPIdata::sout)<<"LocGridOrbitals::add_matrix_elements() for
    // object "
    //             <<mat.name()<<endl;
    fill_dist_matrix_tm_.start();

#ifdef USE_MPI
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sm(comm, dst,
        remote_tasks_DistMatrix_, sparse_distmatrix_nb_tasks_per_partitions_);
#else
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sm(
        0, dst, remote_tasks_DistMatrix_);
#endif

    fillSparseDistMatrix(sm, global_indexes, tol);

    // if(onpe0)(*MPIdata::sout)<<"LocGridOrbitals::add_matrix_elements()"<<endl;
    sm.parallelSumToDistMatrix();

    fill_dist_matrix_tm_.stop();
}

template <class T>
void LocalMatrices<T>::init(const dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
    const vector<vector<int>>& global_indexes,
    const dist_matrix::SubMatricesIndexing<DISTMATDTYPE>& submat_indexing)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();

    dist_matrix::SubMatrices<DISTMATDTYPE> submat(
        "submat", global_indexes, comm, mat, submat_indexing);
    submat.gather(mat);

    const short chromatic_number = (short)global_indexes[0].size();

    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        T* local_mat = ptr_matrices_[iloc];
        for (short icolor = 0; icolor < chromatic_number; icolor++)
            for (short jcolor = 0; jcolor < chromatic_number; jcolor++)
                local_mat[icolor + chromatic_number * jcolor]
                    = (T)submat.val(icolor, jcolor, iloc);
    }
}

template <class T>
void LocalMatrices<T>::applyMask(const LocalMatrices& mask)
{
    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        T* local_mat  = ptr_matrices_[iloc];
        T* local_mask = mask.ptr_matrices_[iloc];
        for (int j = 0; j < n_; j++)
        {
            for (int i = 0; i < m_; i++)
            {
                int index = i + j * m_;
                local_mat[index] *= local_mask[index];
            }
        }
    }
}

template <class T>
void LocalMatrices<T>::setMaskThreshold(
    const T min_threshold, const T max_threshold)
{
    for (short iloc = 0; iloc < subdiv_; iloc++)
    {
        T* local_mat = ptr_matrices_[iloc];
        for (int j = 0; j < n_; j++)
        {
            for (int i = 0; i < m_; i++)
            {
                int index = i + j * m_;
                if (local_mat[index] > max_threshold
                    || local_mat[index] < min_threshold)
                    local_mat[index] = 0.;
                else
                    local_mat[index] = 1.;
            }
        }
    }
}

template class LocalMatrices<double>;
template class LocalMatrices<float>;

template void LocalMatrices<double>::copy(const LocalMatrices<float>& mat);
