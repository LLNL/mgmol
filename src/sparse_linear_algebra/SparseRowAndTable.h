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
#ifndef _SPARSEROWANDTABLE_H_
#define _SPARSEROWANDTABLE_H_

#include "SparseRow.h"
#include "Table.h"

#include <vector>

class SparseRowAndTable : public SparseRow
{

    Table* pos_;

public:
    // constructor
    SparseRowAndTable(const int nnz = DEFAULT_ROW_NNZ) : SparseRow(nnz)
    {
        pos_ = new Table(1);
    }
    // copy constructor
    SparseRowAndTable(const SparseRowAndTable& row) : SparseRow(row)
    {
        pos_ = new Table(1);
        buildRowTable(*pos_);
    }
    SparseRowAndTable(const SparseRow& row) : SparseRow(row)
    {
        pos_ = new Table(1);
        buildRowTable(*pos_);
    }
    // destructor
    ~SparseRowAndTable() override { delete pos_; }
    void assign(std::vector<int>& coldata, std::vector<double>& colvals);
    void assign(const int nnzrow, const int* coldata, const double* colvals);

    void insertEntry(const int col, const double val) override
    {
        //       assert((int *)pos_->get_value(col) == NULL);
        SparseRow::insertEntry(col, val);
        pos_->insert(col);
    }

    int updateRow(
        const int col, const double val, const INSERTMODE mode) override
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

    int updateRowAdd(const int col, const double val) override
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

    int updateRowInsert(const int col, const double val) override
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

    int updateRow(const int count, const int* const cols,
        const double* const vals, const INSERTMODE mode) override;
    int updateRowAdd(const int count, const int* const cols,
        const double* const vals) override;
    int updateRowInsert(const int count, const int* const cols,
        const double* const vals) override;
    void reset();

    /* get value on local row */
    double getColumnEntry(const int col) const override
    {
        int* pos = (int*)pos_->get_value(col);
        if (pos == nullptr)
            return 0.0;
        else
            return SparseRow::getEntryFromPosition(*pos);
    }

    // return column position
    int getColumnPosition(const int col) const override
    {
        int* pos = (int*)pos_->get_value(col);
        return (pos != nullptr) ? *pos : -1;
    }
};

#endif
