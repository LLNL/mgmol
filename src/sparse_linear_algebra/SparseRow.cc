// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SparseRow.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <mpi.h>

int SparseRow::updateRow(const int count, const int* const cols,
    const double* const vals, const INSERTMODE mode)
{
    int newnnz    = 0;
    const int end = (int)vals_.size();
    if (end != 0)
    {
        // sort data
        sortData(0, end - 1);
        /* begin */
        for (int k = 0; k < count; k++)
        {
            // check if column entry exists on local row or not
            int pos = getSortedDataColumnPosition(0, end, cols[k]);
            /* column entry exists. Just insert/ add to it */
            if (pos != -1)
                updateEntry(pos, vals[k], mode);
            else // column entry does not exist. Insert new column entry
            {
                insertEntry(cols[k], vals[k]);
                newnnz++;
            }
        }
    }
    else // row is empty so insert new entry -- equivalent to calling
         // initializeLocalRow()
    {
        assign(count, cols, vals);
        newnnz = count;
    }
    return newnnz;
}

int SparseRow::updateRowAdd(
    const int count, const int* const cols, const double* const vals)
{
    int newnnz    = 0;
    const int end = (int)vals_.size();
    if (end != 0)
    {
        // sort data
        sortData(0, end - 1);
        /* begin */
        for (int k = 0; k < count; k++)
        {
            // check if column entry exists on local row or not
            int pos = getSortedDataColumnPosition(0, end, cols[k]);
            /* column entry exists. Just insert/ add to it */
            if (pos != -1)
                updateEntryAdd(pos, vals[k]);
            else // column entry does not exist. Insert new column entry
            {
                insertEntry(cols[k], vals[k]);
                newnnz++;
            }
        }
    }
    else // row is empty so insert new entry -- equivalent to calling
         // initializeLocalRow()
    {
        assign(count, cols, vals);
        newnnz = count;
    }
    return newnnz;
}

int SparseRow::updateRowInsert(
    const int count, const int* const cols, const double* const vals)
{
    int newnnz    = 0;
    const int end = (int)vals_.size();
    if (end != 0)
    {
        // sort data
        sortData(0, end - 1);
        /* begin */
        for (int k = 0; k < count; k++)
        {
            // check if column entry exists on local row or not
            int pos = getSortedDataColumnPosition(0, end, cols[k]);
            /* column entry exists. Just insert/ add to it */
            if (pos != -1)
                updateEntryInsert(pos, vals[k]);
            else // column entry does not exist. Insert new column entry
            {
                insertEntry(cols[k], vals[k]);
                newnnz++;
            }
        }
    }
    else // row is empty so insert new entry -- equivalent to calling
         // initializeLocalRow()
    {
        assign(count, cols, vals);
        newnnz = count;
    }
    return newnnz;
}

void SparseRow::sortData(const int begin, const int end)
{
    assert((int)cols_.size() > 0);
    assert(begin < (int)cols_.size());
    assert(end < (int)cols_.size());
    //       assert(begin <= end);

    if (end <= begin) return;

    int i   = begin;
    int j   = end;
    int mid = (begin + end) >> 1;
    int x   = ceil(cols_[mid]);

    do
    {
        while (cols_[i] < x)
            i++;
        while (cols_[j] > x)
            j--;
        if (i <= j)
        {
            std::swap(cols_[i], cols_[j]);
            std::swap(vals_[i], vals_[j]);
            i++;
            j--;
        }
    } while (i <= j);
    // recurse
    if (begin < j) sortData(begin, j);
    if (i < end) sortData(i, end);
}

// dot product between row and vector of matrix length
double SparseRow::dotVec(const double* x)
{
    // TODO: Needs an assert on vector length

    double val                          = 0.;
    std::vector<int>::const_iterator it = cols_.begin();
    for (int k = 0; it != cols_.end(); ++it, k++)
    {
        val += vals_[k] * x[*it];
    }

    return val;
}

// compute the pnorm of the row = (sum[ |x_i|^p])^(1/p)
// special case of p=0 returns the max (L-infinity) norm
double SparseRow::pnorm(const int p)
{
    double nrm = 0.;
    if (p > 0)
    {
        for (int k = 0; k < (int)vals_.size(); k++)
        {
            nrm += pow(fabs(vals_[k]), (double)p);
        }
        double xp = 1 / (double)p;
        nrm       = pow(nrm, xp);
    }
    else // return L-infinity norm
    {
        for (int k = 0; k < (int)vals_.size(); k++)
        {
            double aval = fabs(vals_[k]);
            nrm         = nrm > aval ? nrm : aval;
        }
    }
    return nrm;
}
