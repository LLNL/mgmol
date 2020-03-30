// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
 * CSR/CSC matrix used for linear system solver operations
 */

#ifndef _LINEARSOLVERMATRIX_H_
#define _LINEARSOLVERMATRIX_H_

#include "VariableSizeMatrix.h"
#include <vector>

typedef float lsdatatype;

template <class T>
class LinearSolverMatrix
{
    typedef typename std::vector<T>::iterator TvecIterator;
    typedef typename std::vector<T>::const_iterator const_TvecIterator;
    static Timer matrix_insert_tm_;
    static Timer matrix_init_tm_;
    static Timer matrix_initsq_tm_;
    static Timer matvec_tm_;

    int n_; /* matrix size */
    std::vector<int> p_; /* row pointers (of size n+1) */
    std::vector<int> i_; /* column indices, of size nzmax */
    std::vector<T> x_; /* numerical values, of size nzmax */
    std::vector<T> scale_vec_; /* vector of row scaling coeffs */
    bool isrescaled_; /* flag to indicate whether matrix rows have been rescaled
                       */

    void init(
        const int alloc_size, const int nzmax); /* initialize data structs */

public:
    LinearSolverMatrix(const int alloc_size, const int nzmax); // constructor
    void insertMatrixElement(const int pos, const int row, const T val,
        const INSERTMODE mode); /* insert or add matrix entry at position pos */
    void insertMatrixEntry(const int pos, const int row,
        const T val); /* insert/ replace matrix entry at position pos */
    void addToMatrixEntry(
        const int pos, const T val); /* add to matrix entry at position pos */
    void initRow(const int len, const std::vector<int>& cols,
        const std::vector<T>& vals); // initialize a single row of the matrix
    void initRow(const std::vector<int>& cols,
        const std::vector<T>& vals); // initialize a single row of the matrix
    void initRow(const int len, const std::vector<int>::iterator cols_begin,
        const TvecIterator vals_begin); // initialize a single row of the matrix
                                        // - using iterators
    void initRowNNZPattern(const int len,
        const std::vector<int>& cols); // initialize a single row of the matrix
                                       // - nonzero pattern only
    void initRowNNZPattern(
        const std::vector<int>& cols); // initialize a single row of the matrix
                                       // - nonzero pattern only
    void initRowNNZPattern(const int len,
        const std::vector<int>::iterator
            begin); // initialize a single row of the matrix - nonzero pattern
                    // only (with iterators)
    void initializeNonZeroPattern(
        const T val); // Initialize nonzero pattern to take the value val
    void init(const VariableSizeMatrix<sparserow>& vsmat,
        const bool rescale
        = false); // convert from VariableSizeMatrix to LinearSolverMatrix
    void initSquareMat(const VariableSizeMatrix<sparserow>& vsmat,
        const bool rescale = false); // convert from VariableSizeMatrix to
                                     // square LinearSolverMatrix
    void matvec(const std::vector<T>& x,
        std::vector<T>& y) const; // Matrix-vector multiplication operation
    void matvec(const double* const x,
        double* y) const; // Matrix-vector multiplication operation
    void matvec(const float* const x,
        float* y) const; // Matrix-vector multiplication operation
    void matvec(const double* const x,
        float* y) const; // Matrix-vector multiplication operation
    void matvec(const float* const x,
        double* y) const; // Matrix-vector multiplication operation
    T getEntry(const int row, const int col); // get column entry
    ~LinearSolverMatrix(); // destructor

    void reset()
    {
        p_.clear();
        i_.clear();
        x_.clear();
        scale_vec_.clear();
        n_ = 0;
    }

    void reset(const int alloc_size, const int nzmax)
    {
        reset();

        init(alloc_size, nzmax);
    }

    static void printTimers(std::ostream& os)
    {
        matrix_init_tm_.print(os);
        matrix_initsq_tm_.print(os);
        matrix_insert_tm_.print(os);
        matvec_tm_.print(os);
    }

    /* get matrix size */
    int n() const { return n_; }

    /* get nnz */
    int nnz() const { return (int)x_.size(); }

    /* get nonzero pointer value for row i */
    int nzptrval(const int i) const { return p_[i]; }

    /* get column index at position pos */
    int getColumnIndex(const int pos) const { return i_[pos]; }

    /* get value at position pos */
    T getRowEntry(const int pos) const { return x_[pos]; }

    /* Get a vector iterator over non-zero column entries at position pos */
    std::vector<int>::const_iterator getColumnIterator(const int pos) const
    {
        std::vector<int>::const_iterator it = i_.begin() + pos;

        return it;
    }
    /* Get a vector iterator over non-zero entries at position pos */
    TvecIterator getDataIterator(const int pos)
    {
        TvecIterator it = x_.begin() + pos;

        return it;
    }

    /* get pointer to matrix data at position pos */
    T* getPtrToData(const int pos)
    {
        TvecIterator it = getDataIterator(pos);

        return &(*it);
    }

    /* check if matrix is rescaled */
    bool isrescaled() const { return isrescaled_; }

    /* return scaling factor for local row i */
    T getScale(const int i) const { return scale_vec_[i]; }

    /* scale the matrix */
    void scal(const double alpha)
    {
        LinearAlgebraUtils<MemorySpace::Host>::MPscal(n_, alpha, &x_[0]);
    }

    /* ilu factorization - we need templates to allow the precon to be of a
     * different datatype than the matrix */
    template <typename T2>
    void ilu0(LinearSolverMatrix<T2>& L, std::vector<T2>& D,
        LinearSolverMatrix<T2>& U);

    /* Solve with strict lower triangular part of matrix */
    void Lsolve(
        double* const x) const; // In-place triangular solve with L (assume x
                                // contains rhs, and diagonal is all 1)
    void Lsolve(
        float* const x) const; // In-place triangular solve with L (assume x
                               // contains rhs, diagonal is all 1)
    /* Solve with strict upper triangular part of matrix */
    void Usolve(double* const x, const std::vector<T>& diag)
        const; // In-place triangular solve with U (assume x contains rhs,
               // diagonal stored in diag)
    void Usolve(float* const x, const std::vector<T>& diag)
        const; // In-place triangular solve with U (assume x contains rhs,
               // diagonal stored in diag)
};

template <class T>
Timer LinearSolverMatrix<T>::matrix_insert_tm_(
    "LinearSolverMatrix::insertEntry");

template <class T>
Timer LinearSolverMatrix<T>::matrix_init_tm_("LinearSolverMatrix::initialize");

template <class T>
Timer LinearSolverMatrix<T>::matrix_initsq_tm_(
    "LinearSolverMatrix::initialize_sq");

template <class T>
Timer LinearSolverMatrix<T>::matvec_tm_("LinearSolverMatrix::matvec");

#endif
