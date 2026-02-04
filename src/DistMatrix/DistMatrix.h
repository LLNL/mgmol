// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DISTMATRIX_H
#define MGMOL_DISTMATRIX_H

#include "MGmol_scalapack.h"
#include "Timer.h"
#include "mputils.h"

#include <cassert>
#include <complex>
#include <iostream>
#include <string.h>
#include <vector>

#include <mpi.h>

#ifndef MGMOL_DISTMATRIX_ERROR

#define MGMOL_DISTMATRIX_ERROR(X)                                              \
    std::cerr << "ERROR in file " << __FILE__ << " at line " << __LINE__       \
              << std::endl;                                                    \
    std::cerr << "Error Message: " << X << std::endl;                          \
    MPI_Abort(comm_global_, EXIT_FAILURE);

#endif

typedef double DISTMATDTYPE;
// typedef float DISTMATDTYPE;

namespace dist_matrix
{
template <class T>
class DistVector;

class BlacsContext;

template <class T>
class DistMatrix
{
private:
    static Timer allgather_tm_;
    static Timer pdgemr2d_tm_;
    static Timer potri_tm_;
    static Timer potrf_tm_;
    static Timer trtri_tm_;

    static int distmatrix_def_block_size_;
    static BlacsContext* default_bc_;

    std::string object_name_;

    const BlacsContext& bc_;

    const MPI_Comm comm_global_;

    int ictxt_;
    int lld_; // leading dimension of local matrix
    int m_, n_; // size of global matrix
    int mb_, nb_; // size of blocks
    int mloc_, nloc_; // size of local array
    int nprow_, npcol_; // number of process rows and cols in Blacs context
    int myrow_, mycol_; // position of my process in the process grid
    int mblocks_, nblocks_; // number of local blocks
    int desc_[9];
    bool active_;
    bool m_incomplete_, n_incomplete_; // this process has an incomplete block

    void setDiagonalValues(const T* const dmat);

    double dot(const DistMatrix<T>& a) const;

protected:
    // local array
    std::vector<T> val_;

    // size of local data (mloc_ * nloc_)
    int size_;

public:
    std::string name() const { return object_name_; }

    static void setBlockSize(const int nb) { distmatrix_def_block_size_ = nb; }
    static int getBlockSize() { return distmatrix_def_block_size_; }

    static void setDefaultBlacsContext(BlacsContext* bc) { default_bc_ = bc; }

    static void printTimers(std::ostream& os)
    {
        allgather_tm_.print(os);
        pdgemr2d_tm_.print(os);
        potri_tm_.print(os);
        trtri_tm_.print(os);
        potrf_tm_.print(os);
    }

    int nprow() const { return nprow_; }
    int npcol() const { return npcol_; }
    const T val(const int i) const { return val_[i]; }
    void setval(const int i, const T val)
    {
        assert(i < size_);
        val_[i] = val;
    }
    void setval(const int l, const int m, const int x, const int y, const T val)
    {
        val_[l * mb_ + x + (m * nb_ + y) * mloc_] = val;
    }

    // add shift to diagonal, to shift eigenvalues
    void shift(const T shift);

    // input: global index (i,j)
    void addval(const int i, const int j, const T val)
    {
        int dummy    = 0;
        int ig       = i + 1;
        const int il = INDXG2L(&ig, &mb_, &dummy, &dummy, &nprow_) - 1;
        int jg       = j + 1;
        const int jl = INDXG2L(&jg, &nb_, &dummy, &dummy, &npcol_) - 1;
        val_[il + jl * mloc_] += val;
    }

    void setVal(const int i, const int j, const T val);

    int ictxt(void) const { return ictxt_; }
    int lld(void) const { return lld_; }
    int m(void) const { return m_; } // size of global matrix
    int n(void) const { return n_; } // size of global matrix
    int mb(void) const { return mb_; } // size of blocks
    int nb(void) const { return nb_; } // size of blocks
    int mloc(void) const { return mloc_; } // size of local array
    int nloc(void) const { return nloc_; } // size of local array
    int localsize(void) const { return mloc_ * nloc_; } // local size of val
    size_t memsize(void) const { return m_ * n_ * sizeof(T); }
    size_t localmemsize(void) const { return mloc_ * nloc_ * sizeof(T); }
    const int* desc(void) const { return &desc_[0]; }

    // local block size of block (l,m)
    int mbs(const int l) const
    {
        // if l is the last block and this process holds an incomplete
        // block, then the size is m_%mb_. Otherwise, it is mb_.
        return ((l == mblocks_ - 1) && (m_incomplete_)) ? m_ % mb_ : mb_;
    }
    int nbs(const int m) const
    {
        // if m is the last block and this process holds an incomplete
        // block, then the size is n_%nb_. Otherwise, it is nb_.
        return ((m == nblocks_ - 1) && (n_incomplete_)) ? n_ % nb_ : nb_;
    }

    // number of local blocks
    int mblocks(void) const { return mblocks_; }
    int nblocks(void) const { return nblocks_; }

    int myrow(void) const { return myrow_; }
    int mycol(void) const { return mycol_; }

    // functions to compute local indices from global indices

    // index of blocks: element (i,j) is at position (x,y)
    // in local block (l,m) of process (pr,pc)
    int pr(const int i) const
    {
        assert(i < m_);
        assert(i >= 0);
        return (i / mb_) % nprow_;
    }

    int y(const int j) const
    {
        assert(j < n_);
        assert(j >= 0);
        return j % nb_;
    }
    int pc(const int j) const
    {
        assert(j < n_);
        assert(j >= 0);
        return (j / nb_) % npcol_;
    }

    // global indices:
    // (i,j) is the global index of element (x,y) of block (l,m)
    int i(const int l, const int x) const
    {
        return (l * nprow_ + myrow_) * mb_ + x;
    }
    int j(const int m, const int y) const
    {
        return (m * npcol_ + mycol_) * nb_ + y;
    }

    int indxl2grow(const int indxloc) const;
    int indxl2gcol(const int indxloc) const;

    // active flag: the matrix has elements on this process
    bool active(void) const { return active_; }

    DistMatrix<T>(const std::string& name);

    DistMatrix<T>(const std::string& name, const BlacsContext& bc);

    // Construct a DistMatrix of dimensions m,n
    DistMatrix<T>(const std::string& name, const BlacsContext& bc, const int m,
        const int n);

    explicit DistMatrix<T>(const std::string& name, const int m, const int n);

    // Construct a DistMatrix of dimensions m,n
    DistMatrix<T>(const std::string& name, const BlacsContext& bc, const int m,
        const int n, const int mb, const int nb);

    DistMatrix<T>(const std::string& name, const int m, const int n,
        const int mb, const int nb);

    explicit DistMatrix<T>(const std::string& name, const int m);

    // copy constructor
    DistMatrix<T>(const DistMatrix<T>& a);

    // Construct a DistMatrix from a duplicated matrix a, with leading dim lda
    DistMatrix<T>(const std::string& name, const BlacsContext& bc,
        const int lda, const T* const a, const int m, const int n, const int mb,
        const int nb);

    DistMatrix<T>(const std::string& name, const BlacsContext& bc,
        const int lda, const T* const a, const int m, const int n);

    // Construct a diagonal DistMatrix from a vector dmat of diagonal elements
    DistMatrix<T>(const std::string& name, const BlacsContext&,
        const T* const dmat, const int m, const int n, const int mb,
        const int nb);

    // Construct a diagonal DistMatrix from a vector dmat of diagonal elements
    DistMatrix<T>(const std::string& name, const BlacsContext&,
        const T* const dmat, const int m, const int n);
    DistMatrix<T>(const std::string& name, const T* const dmat, const int m);

    DistMatrix<T>& operator=(const DistMatrix<T>& a);
    DistMatrix<T>& assign(const DistMatrix<T>&, const int, const int);
    DistMatrix<T>& operator+=(const DistMatrix<T>& a)
    {
        axpy(1.0, a);
        return *this;
    }
    DistMatrix<T>& operator-=(const DistMatrix<T>& a)
    {
        axpy(-1.0, a);
        return *this;
    }

    ~DistMatrix() { }

    void identity(void);

    // gathers distributed data from all tasks and distribute
    // combined data to all tasks
    void allgather(T* const a, const int lda) const;

    void resize(const int m, const int n, const int mb, const int nb);
    void init(const T* const a, const int lda);
    void add(const T* const a, const int lda);

    void initFromReplicated(double* const a, const int lda);

    void clear(void);

    // compute trace of matrix product
    // trace( (*this)^T*A )
    double traceProduct(const DistMatrix<T>& A) const;

    void axpy(const double alpha, const DistMatrix<T>& a);
    void scal(const double alpha);
    double nrm2() const;
    T asum(DistMatrix<T>& a) const;
    T amax(DistMatrix<T>& a) const;

    double trace(void) const;
    void setDiagonal(const std::vector<T>& diag_values);
    void setRandom(const T minv, const T maxv);

    // matrix * matrix
    // this = alpha*op(A)*op(B)+beta*this
    void gemm(const char transa, const char transb, const T alpha,
        const DistMatrix<T>& a, const DistMatrix<T>& b, const T beta);

    // symmetric_matrix * matrix
    // this = alpha * A * B + beta * this
    void symm(const char side, const char uplo, const T alpha,
        const DistMatrix<T>& a, const DistMatrix<T>& b, const T beta);

    // matrix * vector
    // this = beta*this + alpha * A * x
    void gemv(const char trans, const T alpha, const DistMatrix<T>& a,
        const DistMatrix<T>& x, const T beta);
    void matvec(DistVector<T>& v, DistVector<T>& y);
    // symmetric version of matvec
    void symv(const char uplo, const T alpha, const DistMatrix<T>& a,
        const DistMatrix<T>& x, const T beta);

    // symmetric rank k update
    // this = beta * this + alpha * A * A^T  (trans=='n')
    // this = beta * this + alpha * A^T * A  (trans=='t')
    void syrk(char uplo, char trans, T alpha, DistMatrix<T>& a, T beta);

    // matrix transpose
    // this = alpha * transpose(A) + beta * this
    void transpose(const T alpha, const DistMatrix<T>& a, const T beta);

    void trset(const char);

    void trmm(const char, const char, const char, const char, const T,
        const DistMatrix<T>&);
    void trsm(const char, const char, const char, const char, const T,
        const DistMatrix<T>&);
    void trtrs(const char, const char, const char, DistMatrix<T>&) const;

    // Cholesky decomposition of a symmetric matrix
    int potrf(char uplo);
    // Inverse of a symmetric matrix from Cholesky factor
    int potri(char uplo);
    void potrs(char, DistMatrix<T>&);

    int trtri(char uplo, char diag);

    void getrf(std::vector<int>&);
    void getrs(char, DistMatrix<T>&, std::vector<int>&);
    double norm(char);
    double pocon(char, T);
    void sygst(int, char, const DistMatrix<T>&);
    void syev(char, char, std::vector<T>&, DistMatrix<T>&);
    void sygv(char, char, const char, DistMatrix<T>&, std::vector<T>&,
        DistMatrix<T>&);
    void gesvd(char jobu, char jobvt, std::vector<T>& s, DistMatrix<T>& u,
        DistMatrix<T>& v);

    void getsub(const DistMatrix<T>& a, int m, int n, int ia, int ja);

    int iamax(const int j, T& val);

    void print(std::ostream& os) const;
    void print(
        std::ostream& os, const int, const int, const int, const int) const;
    void printMM(std::ostream& os) const;

    // get value for global index (i,j)
    // assuming we are on the right processor to get it!
    T getVal(const int i, const int j) const;
    T getValOnPE(const int i, const int j) const;

    void assign(const T* const v, const int n)
    {
        if (active_)
        {
            assert(n <= size_);
            memcpy(val_.data(), v, n * sizeof(T));
        }
    }
    void assignColumn(const T* const v, const int i)
    {
        if (active_)
        {
            assert(i <= nloc_);
            if (size_ > 0) memcpy(&val_[i * mloc_], v, mloc_ * sizeof(T));
        }
    }

    // assign data only on first MPI task
    void assign0(const std::vector<T>& v)
    {
        if (myrow_ == 0 && mycol_ == 0)
        {
            assert(v.size() <= static_cast<unsigned int>(size_));
            memcpy(val_.data(), v.data(), v.size() * sizeof(T));
        }
    }

    void swapColumns(const int j1, const int j2);

    void swap(DistMatrix<T>& dm) { val_.swap(dm.val_); }

    void axpyColumn(const int icol, const double alpha, const DistMatrix<T>& x);
    void scalColumn(const int icol, const double alpha);

    void copyDataToVector(std::vector<T>& v) const
    {
        if (active_)
        {
            assert((int)v.size() >= size_);
            memcpy(&v[0], &val_[0], mloc_ * nloc_ * sizeof(T));
        }
    }
    double dotColumns(const int i, const int j) const
    {
        return LinearAlgebraUtils<MemorySpace::Host>::MPdot(
            mloc_, &val_[mloc_ * i], &val_[mloc_ * j]);
    }
    double dotColumns(const int i, const DistMatrix<T>& a, const int j) const
    {
        return LinearAlgebraUtils<MemorySpace::Host>::MPdot(
            mloc_, &val_[mloc_ * i], &a.val_[mloc_ * j]);
    }
    double sumProdElements(const DistMatrix<T>& a) const;
    void getDiagonalValues(T* const dmat) const;

    MPI_Comm comm_global() const;
};

template <class T>
std::ostream& operator<<(std::ostream& os, DistMatrix<T>& a);

template <class T>
int DistMatrix<T>::distmatrix_def_block_size_ = 32;

template <class T>
BlacsContext* DistMatrix<T>::default_bc_ = nullptr;

template <class T>
Timer DistMatrix<T>::allgather_tm_("DistMatrix::matgather");

template <class T>
Timer DistMatrix<T>::pdgemr2d_tm_("DistMatrix::pdgemr2d");

template <class T>
Timer DistMatrix<T>::potri_tm_("DistMatrix::potri");

template <class T>
Timer DistMatrix<T>::potrf_tm_("DistMatrix::potrf");
template <class T>
Timer DistMatrix<T>::trtri_tm_("DistMatrix::trtri");

} // namespace

#endif
