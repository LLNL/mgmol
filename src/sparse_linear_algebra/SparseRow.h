// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
 * Sparse row container for variable size matrix class
 */
#ifndef MGMOL__SPARSEROW_H_
#define MGMOL__SPARSEROW_H_

#include "MGmol_MPI.h"
#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "Table.h"

#include <algorithm>
#include <cassert>
#include <vector>

typedef enum INSERTMODE
{
    INSERT,
    ADD
} INSERTMODE;

/* define maximum local matrix size */
// #define MAX_MAT_SIZE   15000
/* define default tolerance for pruning matrix entries */
#define MAT_TOL 1.0e-14
/* define maximum number of print rows */
#define MAX_PRINT_ROWS 100
/* define default number of print rows for diagnostics */
#define NUM_PRINT_ROWS 5

#define DEFAULT_ROW_NNZ 10

class SparseRow
{
    std::vector<double> vals_;
    std::vector<int> cols_;
    bool issorted_;

    // perform a quicksort of local data. begin and end refer to
    // the lower and upper indexes (inclusive) of the region of the array
    // that is to be sorted
    void sortData(const int begin, const int end);

    std::vector<int>& cols() { return cols_; }
    std::vector<double>& vals() { return vals_; }

public:
    SparseRow(const int nnz = DEFAULT_ROW_NNZ)
    {
        vals_.reserve(nnz);
        cols_.reserve(nnz);
        issorted_ = false;
    }
    SparseRow(const SparseRow& row)
        : vals_(row.vals_), cols_(row.cols_), issorted_(row.issorted_)
    {
    } // copy constructor
    virtual ~SparseRow() { issorted_ = false; } // destructor
    virtual int updateRow(const int count, const int* const cols,
        const double* const vals, const INSERTMODE mode);
    virtual int updateRowAdd(
        const int count, const int* const cols, const double* const vals);
    virtual int updateRowInsert(
        const int count, const int* const cols, const double* const vals);
    /* get value on local row */
    virtual double getColumnEntry(const int col) const
    {
        /*
               if(issorted_)
               {
                  const int end = (int)vals_.size();
                  int pos = getSortedDataColumnPosition(0, end, col);
                  if(pos != -1)
                     return vals_[pos];
                  else
                     return 0.0;
               }
               else
               {
        */
        std::vector<int>::const_iterator it
            = find(cols_.begin(), cols_.end(), col);
        if (it != cols_.end())
        {
            int pos = it - cols_.begin();
            return vals_[pos];
        }
        else
            return 0.0;
        //       }
    }

    void assign(std::vector<int>& coldata, std::vector<double>& colvals)
    {
        const int nnzrow = (int)coldata.size();
        cols_.resize(nnzrow);
        vals_.resize(nnzrow);
        // copy data
        memcpy(&cols_[0], &coldata[0], nnzrow * sizeof(int));
        memcpy(&vals_[0], &colvals[0], nnzrow * sizeof(double));
        issorted_ = false;
    }

    void assign(const int nnzrow, const int* coldata, const double* colvals)
    {
        cols_.resize(nnzrow);
        vals_.resize(nnzrow);
        // copy data
        for (int i = 0; i < nnzrow; i++)
        {
            cols_[i] = coldata[i];
        }
        for (int i = 0; i < nnzrow; i++)
        {
            vals_[i] = colvals[i];
        }
        issorted_ = false;
        //   memcpy(&cols_[0], &coldata[0], nnzrow*sizeof(int));
        //   memcpy(&vals_[0], &colvals[0], nnzrow*sizeof(double));
    }

    void reset()
    {
        vals_.clear();
        cols_.clear();
        issorted_ = false;
    }

    virtual void insertEntry(const int col, const double val)
    {
        assert(find(cols_.begin(), cols_.end(), col) == cols_.end());
        cols_.push_back(col);
        vals_.push_back(val);
        issorted_ = false;
    }

    virtual void updateEntry(const int pos, const double val, INSERTMODE mode)
    {
        assert(pos < (int)vals_.size());

        if (mode == ADD)
            vals_[pos] += val;
        else
            vals_[pos] = val;
    }

    virtual void updateEntryAdd(const int pos, const double val)
    {
        assert(pos < (int)vals_.size());
        vals_[pos] += val;
    }

    virtual void updateEntryInsert(const int pos, const double val)
    {
        assert(pos < (int)vals_.size());
        vals_[pos] = val;
    }

    /* get number of nonzeros for a local row */
    int nnz() const { return (int)vals_.size(); }

    // return column index at position pos
    int getColumnIndex(const int pos) const
    {
        assert(pos < (int)cols_.size());
        return cols_[pos];
    }
    const std::vector<int>& getColumnIndexes() const { return cols_; }
    const std::vector<double>& getColumnEntries() const { return vals_; }
    double* getPtrToColumnEntries() { return &vals_[0]; }
    // return row value at position pos
    double getEntryFromPosition(const int pos) const
    {
        assert(pos >= 0);
        assert(pos < (int)vals_.size());
        return vals_[pos];
    }

    // return column position
    virtual int getColumnPosition(const int col) const
    {
        /*
               if(issorted_)
               {
                  const int end = (int)vals_.size();
                  return getSortedDataColumnPosition(0, end, col);
               }
               else
               {
        */
        std::vector<int>::const_iterator it
            = find(cols_.begin(), cols_.end(), col);
        if (it != cols_.end())
            return (it - cols_.begin());
        else
            return -1;
        //       }
    }

    void copyDataToArray(int* coldata, double* colvals)
    {
        const int nnz = (int)cols_.size();
        memcpy(coldata, &cols_[0], nnz * sizeof(int));
        memcpy(colvals, &vals_[0], nnz * sizeof(double));
    }

    // scale row
    void scale(const double coeff)
    {
        const int len = (int)vals_.size();
        const int one = 1;
        DSCAL(&len, &coeff, &vals_[0], &one);
    }

    // perform daxpy with row data
    void axpy(const int size, const double alpha, double* y)
    {
        assert(size <= (int)vals_.size());
        const int ione = 1;
        DAXPY(&size, &alpha, &vals_[0], &ione, &y[0], &ione);
    }

    void axpy(SparseRow& arow, const double alpha)
    {
        const int m   = arow.nnz();
        const int one = 1;

        std::vector<int> acols(arow.cols());
        std::vector<double> avals(arow.vals());

        DSCAL(&m, &alpha, &avals[0], &one);

        updateRowAdd(m, &acols[0], &avals[0]);
    }

    virtual int updateRow(
        const int col, const double val, const INSERTMODE mode)
    {
        int newentry    = 0;
        const int index = getColumnPosition(col);
        /* column entry exists. Just insert/ add to it */
        if (index != -1)
            updateEntry(index, val, mode);
        else // column entry does not exist. Insert new column entry
        {
            insertEntry(col, val);
            newentry++;
        }

        return newentry;
    }

    virtual int updateRowAdd(const int col, const double val)
    {
        int newentry    = 0;
        const int index = getColumnPosition(col);
        /* column entry exists. Just insert/ add to it */
        if (index != -1)
            updateEntryAdd(index, val);
        else // column entry does not exist. Insert new column entry
        {
            insertEntry(col, val);
            newentry++;
        }

        return newentry;
    }
    virtual int updateRowInsert(const int col, const double val)
    {
        int newentry    = 0;
        const int index = getColumnPosition(col);
        /* column entry exists. Just insert/ add to it */
        if (index != -1)
            updateEntryInsert(index, val);
        else // column entry does not exist. Insert new column entry
        {
            insertEntry(col, val);
            newentry++;
        }

        return newentry;
    }

    // get the column position of column col by searching through column data to
    // see if it exists. Otherwise return -1. begin and end refer to the lower
    // bound and upper bounds of the range of the data to search for col.
    int getSortedDataColumnPosition(
        const int begin, const int end, const int col) const
    {
        assert(static_cast<int>(cols_.size()) > 0);
        assert(begin < static_cast<int>(cols_.size()));
        assert(end <= static_cast<int>(cols_.size()));
        assert(end > 0);
        assert(begin >= 0);

        int pos = -1;
        if (end > begin)
        {
            std::pair<std::vector<int>::const_iterator,
                std::vector<int>::const_iterator>
                bounds;
            bounds = std::equal_range(
                cols_.begin() + begin, cols_.begin() + end, col);
            if (bounds.first != bounds.second)
                pos = bounds.first - cols_.begin();
        }
        return pos;
    }

    int getSortedDataColumnPosition(const int col) const
    {
        int pos       = -1;
        const int end = (int)cols_.size();
        if (end > 0) pos = getSortedDataColumnPosition(0, end, col);

        return pos;
    }

    // perform a quicksort of local data. begin and end refer to
    // the lower and upper indexes (inclusive) of the region of the array
    // that is to be sorted
    void sortData()
    {
        if (issorted_) return;
        const int end = (int)vals_.size();
        if (end > 0)
        {
            sortData(0, end - 1);
            issorted_ = true;
        }
    }

    // insert column indexes into table
    void buildRowTable(Table& pos) { pos.insert(cols_); }

    bool isColSorted() { return issorted_; }

    // dot product between row and vector of matrix length
    double dotVec(const double* x);

    // dot product between row and vector of matrix length
    double dotVec(const std::vector<double>& x) { return dotVec(&x[0]); }

    // compute the pnorm of the row = (sum[ |x_i|^p])^(1/p)
    // special case of p=0 returns the max (L-infinity) norm
    double pnorm(const int p);
};

#endif
