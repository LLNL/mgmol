// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
 * Variable size csr/csc matrix used for data transfer operations
 */
#ifndef MGMOL_VARIABLESIZEMATRIX_H_
#define MGMOL_VARIABLESIZEMATRIX_H_

#include "LocalMatrices.h"
#include "SparseRow.h"
#include "SparseRowAndTable.h"
#include "SquareSubMatrix.h"
#include "Table.h"
#include "VariableSizeMatrixInterface.h"

#include <iostream>
#include <set>
#include <vector>

// typedef enum INSERTMODE  {INSERT, ADD} INSERTMODE;

/* define maximum and minimum local matrix size */
#define MAX_MAT_SIZE 10000
#define MIN_MAT_SIZE 10
/* define default tolerance for pruning matrix entries */
#define MAT_TOL 1.0e-14
/* define maximum number of print rows */
#define MAX_PRINT_ROWS 100
/* define default number of print rows for diagnostics */
#define NUM_PRINT_ROWS 5

class DataDistribution;

/* define matrix row datatype */
typedef SparseRow sparserow;
typedef SparseRowAndTable sparserowtab;

template <class T>
class VariableSizeMatrix : public VariableSizeMatrixInterface
{
    typedef typename std::vector<T*>::iterator TvecIterator;
    typedef typename std::vector<T*>::const_iterator const_TvecIterator;

    const std::string name_;

    int n_; // the dimension of the matrix
    int nzmax_; // max. nnz in each row
    int totnnz_; // total nnz of matrix
    std::vector<int> lvars_; // Local variables in global indices
    Table* table_; // Hash table for holding global, local index pairs

    std::vector<T*> data_;

public:
    VariableSizeMatrix(
        const std::string& name, const int alloc_size); // setup data structures
    VariableSizeMatrix(const VariableSizeMatrix& A,
        const bool copy_table = true); // Copy constructor
    template <class T2>
    VariableSizeMatrix(const VariableSizeMatrix<T2>& A,
        const bool copy_table = true); // Copy constructor

    VariableSizeMatrix<T>& operator=(const VariableSizeMatrix<T>& a);

    /* initialize a local row of the local matrix */
    void updateLocalRowSquareMatrix(const int count, const int lrindex,
        const int* const cols, const double* const vals,
        const INSERTMODE
            mode); /* update current local row by considering only column
                      indices for which there are rows in the local matrix. ie.
                      ensure a square matrix is preserved */
    void insertNewRow(const int ncols, const int row, const int* cols,
        const double* vals,
        const bool append); /* Augment current matrix by inserting a new row */
    /* initialize matrix data from square local matrix object */
    void insertMatrixElements(
        const LocalMatrices<MATDTYPE, MemorySpace::Host>& ss,
        const std::vector<std::vector<int>>& global_indexes, const int numst,
        const double tol = MAT_TOL);
    void insertMatrixElements(
        const SquareSubMatrix<MATDTYPE>& ss, const double tol);
    void sparsify(const std::vector<bool>& keeprow);
    void sparsify(const std::vector<int>& gids);
    void print(std::ostream& os, const std::vector<int>& locfcns,
        int nrows = NUM_PRINT_ROWS) const;
    void printMatCSR(const char* fname); /* print CSR matrix */
    void printMatBlock2(const int gid0, const int gid1, std::ostream& os);
    //    void printMatMM(ofstream& outfile);				/* print
    //    MM matrix */
    void reset(); /* reset CSR matrix to be reused */
    void clear();
    void setupSparseRows(const std::vector<int>&
            rows); /* reset/ initialize matrix with sparse rows */
    void copyData(const VariableSizeMatrix<T>& A,
        const int n); /* Copy data from matrix A. Copies n rows of A */
    void set2Identity(); /* Set matrix to identity */
    ~VariableSizeMatrix() override; // destructor

    void printMat(const char* fname,
        std::vector<int>& lvec); /* print select rows of CSR matrix */

    template <typename T2>
    double AmultSymBdiag(VariableSizeMatrix<T2>* B, const int row);
    double AmultSymB_ij(VariableSizeMatrix<T>* B, const int row,
        const int col); /* compute ij-th entry of A*B */

    double trace(); /* compute the trace of the matrix */
    double trace(const std::vector<int>&
            rows); /* compute the trace of selected rows of the matrix */
    double getTraceDiagProductWithMat(const std::vector<double>&
            ddiag); /* return sum_i ( ddiag[i]*Mat[i][i] ) */
    void copyDataToArray(int* locvars, int* colidx, double* colvals);

    /* get table value */
    void* getTableValue(const int key) const
    {
        return (*table_).get_value(key);
    }
    /* update/ insert key into table */
    void updateTableValue(const int key) { (*table_).insert(key); }
    /* update/ insert key, value into table */
    void updateTableValue(const int key, const int value)
    {
        (*table_).insert(key, value);
    }
    /* get local size */
    int n() const { return n_; }

    /* get nzmax */
    int nzmax() const
    {
        int nzmax = 0;
        const_TvecIterator it;
        for (it = data_.begin(); it != data_.end(); ++it)
            nzmax = nzmax > (int)(*it)->nnz() ? nzmax : (int)(*it)->nnz();

        return nzmax;
    }

    /* get nzmin */
    int nzmin() const
    {
        int nzmin = n_;
        const_TvecIterator it;
        for (it = data_.begin(); it != data_.end(); ++it)
            nzmin = nzmin < (int)(*it)->nnz() ? nzmin : (int)(*it)->nnz();

        return nzmin;
    }

    /* get nzmax of submatrix from row begin to row end */
    int getNzmaxSubmat(const int begin, const int end)
    {
        if (end >= n_) return 0;

        int nzmax = 0;
        for (int i = begin; i <= end; i++)
            nzmax += (int)data_[i]->nnz();

        return nzmax;
    }
    /* get totnnz */
    int nnzmat() const { return totnnz_; }
    /* get number of nonzeros for a local row */
    int nnzrow(const int row) const
    {
        if (row >= n_) return 0;

        return (int)data_[row]->nnz();
    }

    /* get global index of local variable */
    int getLocalVariableGlobalIndex(const int lrindex) const
    {
        return lvars_[lrindex];
    }

    /* get column position on local row. Return -1 if not on local row */
    bool isColumnHere(const int lrindex, const int col) const
    {
        int colpos = data_[lrindex]->getColumnPosition(col);

        if (colpos != -1)
            return true;
        else
            return false;
    }

    /* set pointer to array of global index of local variables */
    int* rowIndexes() { return &lvars_[0]; }

    /* get (global) column index */
    int getColumnIndex(const int lrindex, const int pos) const
    {
        return data_[lrindex]->getColumnIndex(pos);
    }

    void getColumnIndexes(const int lrindex, std::vector<int>& indexes) const
    {
        indexes = data_[lrindex]->getColumnIndexes();
    }

    void getAllColumnIndexes(std::vector<int>& indexes) const;

    /* get value on local row */
    double getRowEntry(const int lrindex, const int pos) const
    {
        assert(lrindex < n_);
        return data_[lrindex]->getEntryFromPosition(pos);
    }

    void getRowEntries(const int lrindex, std::vector<double>& values) const
    {
        assert(lrindex < n_);

        values = data_[lrindex]->getColumnEntries();
    }

    int getMaxAbsOffDiagonalRowEntry(const int gid, double& value) const;

    int getColumnPos(const int lrindex, const int col)
    {
        return data_[lrindex]->getColumnPosition(col);
    }

    void row_daxpy(
        const int lrindex, const int size, const double alpha, double* y)
    {
        data_[lrindex]->axpy(size, alpha, y);
    }

    Table* getTable() { return table_; }

    /* initialize a local row of the local matrix */
    /* Assumes nnzrow is initially zero - matrix has been reset */
    void initializeLocalRow(
        const int ncols, const int lrindex, const int* cols, const double* vals)
    {
        if (ncols)
        {
            data_[lrindex]->assign(ncols, cols, vals);
            /* update local matrix variables */
#ifdef _OPENMP
#pragma omp atomic
#endif
            totnnz_ += ncols;
        }

        return;
    }

    /* Update current local rows by adding or inserting new columns. */
    void updateLocalRow(const int count, const int lrindex,
        const int* const cols, const double* const vals, const INSERTMODE mode)
    {
        // updateRow_tm_.start();
        totnnz_ += data_[lrindex]->updateRow(count, cols, vals, mode);
        // updateRow_tm_.stop();
        return;
    }
    void updateLocalRowAdd(const int count, const int lrindex,
        const int* const cols, const double* const vals)
    {
        // updateRow_tm_.start();

        const int newnnz = data_[lrindex]->updateRowAdd(count, cols, vals);
#ifdef _OPENMP
#pragma omp atomic
#endif
        totnnz_ += newnnz;
        // updateRow_tm_.stop();

        return;
    }
    void updateLocalRowInsert(const int count, const int lrindex,
        const int* const cols, const double* const vals)
    {
        // updateRow_tm_.start();
        totnnz_ += data_[lrindex]->updateRowInsert(count, cols, vals);
        // updateRow_tm_.stop();
        return;
    }
    /* Update current local row by adding or inserting a new column. */
    void updateLocalRow(const int lrindex, const int col, const double val,
        const INSERTMODE mode)
    {
        // updateRow_tm_.start();
        totnnz_ += data_[lrindex]->updateRow(col, val, mode);
        // updateRow_tm_.stop();
        return;
    }
    void updateLocalRowAdd(const int lrindex, const int col, const double val)
    {
        // updateRow_tm_.start();
        totnnz_ += data_[lrindex]->updateRowAdd(col, val);
        // updateRow_tm_.stop();
        return;
    }
    void updateLocalRowInsert(
        const int lrindex, const int col, const double val)
    {
        // updateRow_tm_.start();
        totnnz_ += data_[lrindex]->updateRowInsert(col, val);
        // updateRow_tm_.stop();
        return;
    }

    /* Update current local entry by adding or inserting a new value.
     * Assumes that the local row index and column position of the entry
     * is known.
     */
    void updateLocalEntry(const int lrindex, const int pos, const double val,
        const INSERTMODE mode)
    {
        /* begin */
        /* Add or insert entry */
        data_[lrindex]->updateEntry(pos, val, mode);

        return;
    }

    void updateLocalEntryAdd(const int lrindex, const int pos, const double val)
    {
        /* begin */
        /* Add or insert entry */
        data_[lrindex]->updateEntryAdd(pos, val);

        return;
    }

    void updateLocalEntryInsert(
        const int lrindex, const int pos, const double val)
    {
        /* begin */
        /* Add or insert entry */
        data_[lrindex]->updateEntryInsert(pos, val);

        return;
    }

    /* Insert entry into matrix */
    void insertMatrixElement(const int row, const int col, const double val,
        const INSERTMODE mode, const bool append)
    {
#ifdef _OPENMP
#pragma omp critical(insertMatrixElement)
#endif
        {
            /* begin */
            /* check if row exists */
            int* rindex = (int*)getTableValue(row);

            if (rindex != nullptr) /* row exists */
            {
                /* insert column */
                updateLocalRow(*rindex, col, val, mode);
            }
            else /* insert new row */
            {
                insertNewRow(1, row, &col, &val, append);
            }
        }
    }

    /* get matrix entry */
    double get_value(const int row, const int col) const
    {
        double value = 0.0;

        int* rindex = (int*)getTableValue(row);
        if (rindex != nullptr) value = data_[*rindex]->getColumnEntry(col);

        return value;
    }
    /* get matrix entries from a local row = lrindex */
    void getLocalRowValues(const int lrindex, const std::vector<int>& cols,
        std::vector<double>& vals) const
    {
        vals.reserve(cols.size());
        for (std::vector<int>::const_iterator it = cols.begin();
             it != cols.end(); ++it)
        {
            vals.push_back(data_[lrindex]->getColumnEntry(*it));
        }
        assert(vals.size() == cols.size());
    }

    /* get matrix entries from a global row*/
    void getRowValues(const int row, const std::vector<int>& cols,
        std::vector<double>& vals) const
    {
        int* rindex = (int*)getTableValue(row);
        if (rindex != nullptr)
        {
            const int lrindex = *rindex;
            getLocalRowValues(lrindex, cols, vals);
            /*
                   T* data=data_[*rindex];
                   vals.reserve(cols.size());

                   for(std::vector<int>::const_iterator it =cols.begin();
                                                        it!=cols.end();
                                                      ++it)
                   {
                       vals.push_back( data->getColumnEntry(*it) );
                   }
            */
        }
        else
        {
            vals.resize(cols.size());
            memset(&vals[0], 0, vals.size() * sizeof(double));
        }

        assert(vals.size() == cols.size());
    }

    /* get matrix entries from a sorted row*/
    void getSortedRowValues(const int row, const std::vector<int>& cols,
        std::vector<double>& vals) const
    {
        int* rindex = (int*)getTableValue(row);
        if (rindex != nullptr)
        {
            T* data = data_[*rindex];
            vals.reserve(cols.size());

            sort_col_tm_.start();
            data_[*rindex]->sortData();
            sort_col_tm_.stop();

            for (std::vector<int>::const_iterator it = cols.begin();
                 it != cols.end(); ++it)
            {
                int pos = data->getSortedDataColumnPosition(*it);
                if (pos != -1)
                    vals.push_back(data->getEntryFromPosition(pos));
                else
                    vals.push_back(0.0);
            }
        }
        else
        {
            vals.resize(cols.size());
            memset(&vals[0], 0, vals.size() * sizeof(double));
        }

        assert(vals.size() == cols.size());
    }

    /* Scale the row of the CSR matrix */
    void scaleRow(const int row, const double coeff)
    {
        int* rindex = (int*)getTableValue(row);
        if (rindex == nullptr) return;

        data_[*rindex]->scale(coeff);
    }

    /* Scale the CSR matrix */
    void scale(const double coeff)
    {

        const int n = n_;
        for (int lrindex = 0; lrindex < n; lrindex++)
            data_[lrindex]->scale(coeff);
    }

    // matrix multiplication operations (locally centered contributions only)
    // flag== true => compute entries for specific nonzero pattern only
    void AmultSymBLocal(VariableSizeMatrix<T>* B, VariableSizeMatrix<T>& C,
        const std::vector<int>& locfcns,
        VariableSizeMatrix<SparseRowAndTable>& pattern, bool flag = true);

    // matrix multiplication operations
    void AmultSymB(VariableSizeMatrix<T>* B, VariableSizeMatrix<T>& C,
        VariableSizeMatrix<SparseRowAndTable>& pattern, bool flag = true);

    const std::vector<int>& lvars() const { return lvars_; }

    // get reference to local row at index rindex
    T& getRow(int rindex) const { return *data_[rindex]; }

    void sortColumnIndexes()
    {
        sort_col_tm_.start();

        for (const_TvecIterator it = data_.begin(); it != data_.end(); ++it)
            (*it)->sortData();
        sort_col_tm_.stop();
    }

    // get pointer to row data
    double* getRowEntries(const int lrindex)
    {
        assert(lrindex < n_);

        return data_[lrindex]->getPtrToColumnEntries();
    }

    void axpy(const double alpha, const VariableSizeMatrix<T>& B);
    void gemv(const double alpha, const std::vector<double>& x,
        const double beta, std::vector<double>& y);

    // compute dot product of matrix row with an array
    double rowDotVec(const int row, const double* x)
    {
        return data_[row]->dotVec(x);
    }

    double pnorm(const int row, const int p) { return data_[row]->pnorm(p); }

    VariableSizeMatrix<T>& operator+=(const VariableSizeMatrix<T>& a)
    {
        axpy(1.0, a);
        return *this;
    }
    VariableSizeMatrix<T>& operator-=(const VariableSizeMatrix<T>& a)
    {
        axpy(-1.0, a);
        return *this;
    }
    std::string name() { return name_; }
};

#endif
