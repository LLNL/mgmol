// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SparseRowAndTable.h"
#include "Table.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>

void SparseRowAndTable::assign(
    std::vector<int>& coldata, std::vector<double>& colvals)
{
    // copy data
    SparseRow::assign(coldata, colvals);
    // setup table for position data
    pos_->insert(coldata);
}

void SparseRowAndTable::assign(
    const int nnzrow, const int* coldata, const double* colvals)
{
    // copy data
    SparseRow::assign(nnzrow, coldata, colvals);
    // setup table for position data
    buildRowTable(*pos_);
}

int SparseRowAndTable::updateRow(const int count, const int* const cols,
    const double* const vals, const INSERTMODE mode)
{
    int newnnz    = 0;
    const int end = nnz();
    if (end != 0)
    {
        /* begin */
        for (int k = 0; k < count; k++)
        {
            // check if column entry exists on local row or not
            int pos = getColumnPosition(cols[k]);
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

int SparseRowAndTable::updateRowAdd(
    const int count, const int* const cols, const double* const vals)
{
    int newnnz    = 0;
    const int end = nnz();
    if (end != 0)
    {
        /* begin */
        for (int k = 0; k < count; k++)
        {
            // check if column entry exists on local row or not
            int pos = getColumnPosition(cols[k]);
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

int SparseRowAndTable::updateRowInsert(
    const int count, const int* const cols, const double* const vals)
{
    int newnnz    = 0;
    const int end = nnz();
    if (end != 0)
    {
        /* begin */
        for (int k = 0; k < count; k++)
        {
            // check if column entry exists on local row or not
            int pos = getColumnPosition(cols[k]);
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

void SparseRowAndTable::reset()
{
    SparseRow::reset();
    pos_->reset();
}
