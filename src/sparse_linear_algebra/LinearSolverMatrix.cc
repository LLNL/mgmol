// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cassert>
#include <vector>

#include "LinearSolverMatrix.h"
#include "VariableSizeMatrix.h"

template <class T>
LinearSolverMatrix<T>::LinearSolverMatrix(const int n, const int nzmax)
{
    init(n, nzmax);
}
template <class T>
LinearSolverMatrix<T>::~LinearSolverMatrix()
{
}

/* insert matrix entry at position pos.
 * Note: this function assumes that the entry
 * already exists and is being changed
 */
template <class T>
void LinearSolverMatrix<T>::insertMatrixElement(
    const int pos, const int row, const T val, const INSERTMODE mode)
{
    matrix_insert_tm_.start();
    assert(pos < (int)x_.size());

    i_[pos] = row;
    if (mode == ADD)
        x_[pos] += val;
    else
        x_[pos] = val;

    matrix_insert_tm_.stop();
    return;
}

/* insert matrix entry at position pos.
 * Note: this function assumes that the entry
 * already exists and is being changed
 */
template <class T>
void LinearSolverMatrix<T>::insertMatrixEntry(
    const int pos, const int row, const T val)
{
    matrix_insert_tm_.start();
    assert(pos < (int)x_.size());

    i_[pos] = row;
    x_[pos] = val;

    matrix_insert_tm_.stop();
    return;
}

/* Add matrix entry at position pos.
 * Note: this function assumes that the entry
 * already exists and is being modified
 */
template <class T>
void LinearSolverMatrix<T>::addToMatrixEntry(const int pos, const T val)
{
    matrix_insert_tm_.start();
    assert(pos < (int)x_.size());

    x_[pos] += val;

    matrix_insert_tm_.stop();
    return;
}

/* initialize a single row of the matrix  */
template <class T>
void LinearSolverMatrix<T>::initRow(
    const int len, const std::vector<int>& cols, const std::vector<T>& vals)
{
    /* begin */
    matrix_insert_tm_.start();
    /* update nonzero column pointer array */
    p_.push_back(p_.back() + len);

    if (len > 0)
    {
        i_.insert(i_.end(), cols.begin(), (cols.begin() + len));
        x_.insert(x_.end(), vals.begin(), (vals.begin() + len));
    }
    matrix_insert_tm_.stop();
    return;
}

/* initialize a single row of the matrix  */
template <class T>
void LinearSolverMatrix<T>::initRow(
    const std::vector<int>& cols, const std::vector<T>& vals)
{
    /* begin */
    matrix_insert_tm_.start();
    const int len = (int)cols.size();
    initRow(len, cols, vals);
    matrix_insert_tm_.stop();
    return;
}

/* initialize a single row of the matrix - (iterator version) */
template <class T>
void LinearSolverMatrix<T>::initRow(const int len,
    const std::vector<int>::iterator cols_begin, const TvecIterator vals_begin)
{
    /* begin */
    matrix_insert_tm_.start();
    /* update nonzero column pointer array */
    p_.push_back(p_.back() + len);

    if (len > 0)
    {
        i_.insert(i_.end(), cols_begin, (cols_begin + len));
        x_.insert(x_.end(), vals_begin, (vals_begin + len));
    }
    matrix_insert_tm_.stop();
    return;
}

/* initialize a single row of the matrix - Non-zero pattern only */
template <class T>
void LinearSolverMatrix<T>::initRowNNZPattern(
    const int len, const std::vector<int>& cols)
{
    /* begin */
    matrix_insert_tm_.start();
    /* update nonzero column pointer array */
    p_.push_back(p_.back() + len);

    if (len > 0)
    {
        i_.insert(i_.end(), cols.begin(), (cols.begin() + len));
        //       x_.insert(x_.end(), vals.begin(), (vals.begin()+len));
    }
    matrix_insert_tm_.stop();
    return;
}

/* initialize a single row of the matrix - Non-zero pattern only */
template <class T>
void LinearSolverMatrix<T>::initRowNNZPattern(const std::vector<int>& cols)
{
    /* begin */
    matrix_insert_tm_.start();
    const int len = (int)cols.size();
    initRowNNZPattern(len, cols);
    matrix_insert_tm_.stop();
    return;
}

/* initialize a single row of the matrix - Non-zero pattern only (iterator
 * version) */
template <class T>
void LinearSolverMatrix<T>::initRowNNZPattern(
    const int len, const std::vector<int>::iterator begin)
{
    /* begin */
    matrix_insert_tm_.start();
    /* update nonzero column pointer array */
    p_.push_back(p_.back() + len);

    if (len > 0)
    {
        i_.insert(i_.end(), begin, (begin + len));
        //       x_.insert(x_.end(), vals.begin(), (vals.begin()+len));
    }
    matrix_insert_tm_.stop();
    return;
}

template <class T>
void LinearSolverMatrix<T>::initializeNonZeroPattern(const T val)
{
    matrix_insert_tm_.start();
    const int nzmax = (int)i_.size();
    x_.resize(nzmax, val);
    matrix_insert_tm_.stop();
    return;
}

/* initialize matrix with data from VariableSizeMatrix object */
template <class T>
void LinearSolverMatrix<T>::init(
    const VariableSizeMatrix<sparserow>& vsmat, const bool rescale)
{
    const int n   = n_;
    const T dfact = 1.0e-2;
    /* begin */
    matrix_init_tm_.start();
    if (rescale)
    {
        T dnew        = 1.0e-2;
        const int one = 1;
        for (int i = 0; i < n; i++)
        {
            const int nnzrow = vsmat.nnzrow(i);
            int len          = nnzrow;
            bool dflag       = false;
            for (int j = 0; j < nnzrow; j++)
            {
                const int key = vsmat.getColumnIndex(i, j);
                int* cindex   = (int*)vsmat.getTableValue(key);
                if (cindex == nullptr) continue;
                const int lcindex = *cindex;
                i_.push_back(lcindex);
                T val = (T)vsmat.getRowEntry(i, j);
                if (lcindex == i)
                {
                    dflag = true;
                    // dfact = val;
                    if (fabs(val) < dfact)
                    {
                        T scale = 1.0 / val;
                        scale *= dnew;
                        // val *= dnew;
                        scale_vec_.push_back(scale);
                    }
                    else
                    {
                        T scale = 1.0; // 1.0/val;
                        scale_vec_.push_back(scale);
                    }
                }
                x_.push_back(val);
            }
            if (dflag == false) // rescale zero diagonals
            {
                i_.push_back(i);
                x_.push_back(dnew);
                scale_vec_.push_back(1.0);
                len = len + 1;
            }
            const int k = (int)x_.size();
            p_.push_back(k);

            /* rescale row */
            DSCAL(&len, (double*)&scale_vec_[i], (double*)getPtrToData(p_[i]),
                &one);
        }
    }
    else
    {
        x_.reserve(n);
        i_.reserve(n);
        p_.reserve(n);
        for (int i = 0; i < n; i++)
        {
            const int nnzrow = vsmat.nnzrow(i);
            bool dflag       = false;
            for (int j = 0; j < nnzrow; j++)
            {
                const int key = vsmat.getColumnIndex(i, j);
                int* cindex   = (int*)vsmat.getTableValue(key);
                if (cindex == nullptr) continue;
                const int lcindex = *cindex;
                i_.push_back(lcindex);
                T val = (T)vsmat.getRowEntry(i, j);
                if (lcindex == i) // rescale small diagonals
                {
                    dflag = true;
                    val   = fabs(val) > dfact ? val : dfact;
                }
                x_.push_back(val);
            }
            if (dflag == false) // rescale zero diagonals
            {
                i_.push_back(i);
                x_.push_back(dfact);
            }
            p_.push_back((int)x_.size());
        }
    }
    matrix_init_tm_.stop();

    isrescaled_ = rescale;
    return;
}

/* initialize square (sub)matrix with data from VariableSizeMatrix object */
template <class T>
void LinearSolverMatrix<T>::initSquareMat(
    const VariableSizeMatrix<sparserow>& vsmat, const bool rescale)
{
    const int n = n_;
    T dfact     = 1.0e-2;
    /* begin */
    matrix_initsq_tm_.start();
    if (rescale == true)
    {
        T dnew        = 1.0e-2;
        const int one = 1;
        for (int i = 0; i < n; i++)
        {
            const int nnzrow = vsmat.nnzrow(i);
            int len          = nnzrow;
            bool dflag       = false;
            for (int j = 0; j < nnzrow; j++)
            {
                const int key = vsmat.getColumnIndex(i, j);
                int* cindex   = (int*)vsmat.getTableValue(key);
                if (cindex == nullptr) continue;
                const int lcindex = *cindex;
                if (lcindex >= n) continue;
                i_.push_back(lcindex);
                T val = (T)vsmat.getRowEntry(i, j);
                if (lcindex == i)
                {
                    dflag = true;
                    // dfact = val;
                    if (fabs(val) < dfact)
                    {
                        T scale = 1.0 / val;
                        scale *= dnew;
                        // val *= dnew;
                        scale_vec_.push_back(scale);
                    }
                    else
                    {
                        T scale = 1.0; // 1.0/val;
                        scale_vec_.push_back(scale);
                    }
                }
                x_.push_back(val);
            }
            if (dflag == false) // rescale zero diagonals
            {
                i_.push_back(i);
                x_.push_back(dnew);
                scale_vec_.push_back(1.0);
                len = len + 1;
            }
            const int k = (int)x_.size();
            p_.push_back(k);

            /* rescale row */
            DSCAL(&len, (double*)&scale_vec_[i], (double*)getPtrToData(p_[i]),
                &one);
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            const int nnzrow = vsmat.nnzrow(i);
            bool dflag       = false;
            for (int j = 0; j < nnzrow; j++)
            {
                const int key = vsmat.getColumnIndex(i, j);
                int* cindex   = (int*)vsmat.getTableValue(key);
                if (cindex == nullptr) continue;
                const int lcindex = *cindex;
                if (lcindex >= n) continue;
                i_.push_back(lcindex);
                T val = (T)vsmat.getRowEntry(i, j);
                if (lcindex == i) // rescale small diagonals
                {
                    dflag = true;
                    val   = fabs(val) > dfact ? val : dfact;
                }
                x_.push_back(val);
            }
            if (dflag == false) // rescale zero diagonals
            {
                i_.push_back(i);
                x_.push_back(dfact);
            }
            const int k = (int)x_.size();
            p_.push_back(k);
        }
    }
    matrix_initsq_tm_.stop();

    isrescaled_ = rescale;
    return;
}

/* matrix vector multiply */
template <class T>
void LinearSolverMatrix<T>::matvec(
    const std::vector<T>& v, std::vector<T>& w) const
{
    matvec_tm_.start();

    const int siz = (int)v.size();
    std::fill_n(w.begin(), siz, 0.);
    //    w.clear();
    //    w.resize(siz);

    for (int i = 0; i < siz; i++)
    {
        const int start = p_[i];
        const int stop  = p_[i + 1];
        for (int k = start; k < stop; k++)
        {
            w[i_[k]] += x_[k] * v[i];
        }
    }
    matvec_tm_.stop();
    return;
}

/* matrix vector multiply */
template <class T>
void LinearSolverMatrix<T>::matvec(const double* const v, double* w) const
{
    matvec_tm_.start();
    const int n = n_;

    memset(w, 0, n * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        const int start                      = p_[i];
        const int stop                       = p_[i + 1];
        std::vector<int>::const_iterator it  = getColumnIterator(start);
        std::vector<int>::const_iterator end = getColumnIterator(stop);
        const double fact                    = v[i];

        for (int k = start; it != end; ++it, k++)
        {
            w[*it] += x_[k] * fact;
        }
    }
    matvec_tm_.stop();
    return;
}

/* matrix vector multiply */
template <class T>
void LinearSolverMatrix<T>::matvec(const float* const v, float* w) const
{
    matvec_tm_.start();
    const int n = n_;

    // for accumulation in double
    double* z = new double[n];
    memset(z, 0, n * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        const int start                      = p_[i];
        const int stop                       = p_[i + 1];
        std::vector<int>::const_iterator it  = getColumnIterator(start);
        std::vector<int>::const_iterator end = getColumnIterator(stop);
        const double fact                    = (double)v[i];

        for (int k = start; it != end; ++it, k++)
        {
            z[*it] += ((double)x_[k] * fact);
        }
    }

    for (int i = 0; i < n; i++)
        w[i] = (float)z[i];

    delete[] z;
    matvec_tm_.stop();
    return;
}

/* matrix vector multiply */
template <class T>
void LinearSolverMatrix<T>::matvec(const double* const v, float* w) const
{
    matvec_tm_.start();
    const int n = n_;

    memset(w, 0, n * sizeof(float));

    for (int i = 0; i < n; i++)
    {
        const int start                      = p_[i];
        const int stop                       = p_[i + 1];
        std::vector<int>::const_iterator it  = getColumnIterator(start);
        std::vector<int>::const_iterator end = getColumnIterator(stop);
        const double fact                    = v[i];

        for (int k = start; it != end; ++it, k++)
        {
            w[*it] += (float)(x_[k] * fact);
        }
    }
    matvec_tm_.stop();
    return;
}

/* matrix vector multiply */
template <class T>
void LinearSolverMatrix<T>::matvec(const float* const v, double* w) const
{
    matvec_tm_.start();
    const int n = n_;

    memset(w, 0, n * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        const int start                      = p_[i];
        const int stop                       = p_[i + 1];
        std::vector<int>::const_iterator it  = getColumnIterator(start);
        std::vector<int>::const_iterator end = getColumnIterator(stop);
        const double fact                    = (double)v[i];

        for (int k = start; it != end; ++it, k++)
        {
            w[*it] += ((double)x_[k] * fact);
        }
    }
    matvec_tm_.stop();
    return;
}

/* get column entry */
template <class T>
T LinearSolverMatrix<T>::getEntry(const int row, const int col)
{
    std::vector<int>::const_iterator start = i_.begin() + p_[row];
    std::vector<int>::const_iterator stop  = i_.begin() + p_[row + 1];

    std::vector<int>::const_iterator it = find(start, stop, col);
    if (it != stop)
    {
        int pos = it - i_.begin();
        return x_[pos];
    }
    else
        return 0.0;
}

template <class T>
void LinearSolverMatrix<T>::Lsolve(double* const x) const
{
    for (int i = 0; i < n_; i++)
    {
        const int k1                         = p_[i];
        const int k2                         = p_[i + 1];
        std::vector<int>::const_iterator row = i_.begin() + k1;
        std::vector<int>::const_iterator end = i_.begin() + k2;
        const_TvecIterator dptr              = x_.begin() + k1;

        const double fact = x[i];
        for (; row != end; ++row, ++dptr)
            x[*row] -= (*dptr) * fact;
    }

    return;
}

template <class T>
void LinearSolverMatrix<T>::Lsolve(float* const x) const
{
    double* z = new double[n_];
    for (int i = 0; i < n_; i++)
        z[i] = (double)x[i];

    for (int i = 0; i < n_; i++)
    {
        const int k1                         = p_[i];
        const int k2                         = p_[i + 1];
        std::vector<int>::const_iterator row = i_.begin() + k1;
        std::vector<int>::const_iterator end = i_.begin() + k2;
        const_TvecIterator dptr              = x_.begin() + k1;

        const double fact = z[i];
        for (; row != end; ++row, ++dptr)
            z[*row] -= ((double)(*dptr) * fact);
    }
    for (int i = 0; i < n_; i++)
        x[i] = (float)z[i];
    delete[] z;
    return;
}

template <class T>
void LinearSolverMatrix<T>::Usolve(
    double* const x, const std::vector<T>& diag) const
{
    assert((int)diag.size() == n_);
    for (int i = n_ - 1; i >= 0; i--)
    {
        x[i] *= diag[i];
        const int k1                         = p_[i];
        const int k2                         = p_[i + 1];
        std::vector<int>::const_iterator row = i_.begin() + k1;
        std::vector<int>::const_iterator end = i_.begin() + k2;
        const_TvecIterator dptr              = x_.begin() + k1;

        const double fact = x[i];
        for (; row != end; ++row, ++dptr)
            x[*row] -= (*dptr) * fact;
    }

    return;
}

template <class T>
void LinearSolverMatrix<T>::Usolve(
    float* const x, const std::vector<T>& diag) const
{
    assert((int)diag.size() == n_);

    double* z = new double[n_];
    for (int i = 0; i < n_; i++)
        z[i] = (double)x[i];

    for (int i = n_ - 1; i >= 0; i--)
    {
        z[i] *= (double)diag[i];
        const int k1                         = p_[i];
        const int k2                         = p_[i + 1];
        std::vector<int>::const_iterator row = i_.begin() + k1;
        std::vector<int>::const_iterator end = i_.begin() + k2;
        const_TvecIterator dptr              = x_.begin() + k1;

        const double fact = z[i];
        for (; row != end; ++row, ++dptr)
            z[*row] -= ((double)(*dptr) * fact);
    }
    for (int i = 0; i < n_; i++)
        x[i] = (float)z[i];
    delete[] z;
    return;
}

template <class T>
void LinearSolverMatrix<T>::init(const int alloc_size, const int nzmax)
{
    n_ = alloc_size;
    if (n_ > 0)
    {
        p_.reserve((n_ + 1));
        p_.push_back(0); /* initialize first entry of nonzero pointer */

        i_.reserve(nzmax);
        x_.reserve(nzmax);
    }
}

template class LinearSolverMatrix<float>;
template class LinearSolverMatrix<double>;
