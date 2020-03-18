// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "PackedCommunicationBuffer.h"
#include "MPIdata.h"
#include "VariableSizeMatrix.h"

#include "../Control.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <mpi.h>

Timer PackedCommunicationBuffer::pack_local_data_tm_(
    "PackedCommunicationBuffer::packLocalData");
Timer PackedCommunicationBuffer::merge_or_insert_new_rows_tm_(
    "PackedCommunicationBuffer::mergeOrInsertNewRows");
Timer PackedCommunicationBuffer::merge_to_existing_rows_tm_(
    "PackedCommunicationBuffer::mergeToExistingRows");
Timer PackedCommunicationBuffer::copy_and_insert_new_rows_tm_(
    "PackedCommunicationBuffer::copyAndInsertNewRows");

int PackedCommunicationBuffer::storage_size_  = 0;
char* PackedCommunicationBuffer::storage_     = nullptr;
char* PackedCommunicationBuffer::send_buffer_ = nullptr;
char* PackedCommunicationBuffer::recv_buffer_ = nullptr;

PackedCommunicationBuffer::PackedCommunicationBuffer(const int bsize)
{
    Control& ct = *(Control::instance());

    if (bsize * 2 > storage_size_)
    {
        if (onpe0 && ct.verbose > 0)
            std::cout << "PackedCommunicationBuffer: resize buffers to 2 x "
                      << bsize << std::endl;
        storage_size_ = 2 * bsize;

        if (storage_ != nullptr) delete[] storage_;
        storage_ = new char[storage_size_];

        send_buffer_ = &storage_[0];
        recv_buffer_ = &storage_[bsize];
    }
    memset(storage_, 0, storage_size_ * sizeof(char));

    packed_nnzrow_ptr_ = nullptr;
    packed_lvars_ptr_  = nullptr;
    packed_pj_ptr_     = nullptr;
    packed_pa_ptr_     = nullptr;
    packed_data_size_  = 0;
}

/* Set data pointer positions on recv buffer */
void PackedCommunicationBuffer::setupPackedDataPointers(const char* data)
{
    assert(data != nullptr);

    // setup_data_ptr_tm_.start();
    int* iptr        = (int*)data;
    packed_num_rows_ = *(iptr++);

    /* check if matrix is nonzero */
    if (packed_num_rows_ != 0)
    {
        const int packed_double_pos = *(iptr++);
        packed_nnzrow_ptr_          = iptr;
        packed_lvars_ptr_ = packed_nnzrow_ptr_ + (packed_num_rows_ + 1);
        packed_pj_ptr_    = packed_lvars_ptr_ + packed_num_rows_;
        packed_pa_ptr_    = (double*)(data + packed_double_pos);

        /* compute size of data in buffer */
        packed_data_size_
            = packed_double_pos
              + packed_nnzrow_ptr_[packed_num_rows_] * sizeof(double);
    }
    else
    {
        /* The size is equal to the size of an int if nrows == 0. */
        packed_data_size_ = sizeof(int);
    }
    // setup_data_ptr_tm_.stop();
}

/* Pack local data into send buffer */
void PackedCommunicationBuffer::packLocalData(
    VariableSizeMatrix<sparserow>& lmat, const int* pos)
{
    assert(lmat.n() * sizeof(double) < (storage_size_ / 2) * sizeof(char));

    pack_local_data_tm_.start();

    /* Get matrix dimension */
    const int lsize = lmat.n();

    /* Manually pack buffer */
    /* check for zero data - base case */
    if (lsize == 0)
    {
        int* const iptr = (int*)&send_buffer_[pos[0]];
        *(iptr)         = 0;
    }
    else
    {
        /* pack matrix size - n */
        int* iptr = (int*)&send_buffer_[pos[0]];
        *(iptr++) = lsize;
        /* pack start position of doubles (matrix entries) */
        *(iptr++) = (int)pos[1];
        /* pack nb. nonzero for each row */
        *(iptr++) = 0;
        int k     = 0;
        for (int i = 0; i < lsize; i++)
        {
            k += lmat.nnzrow(i);
            *(iptr++) = k;
        }
        /* pack matrix data */
        double* dptr = (double*)&send_buffer_[pos[1]];
        lmat.copyDataToArray(iptr, iptr + lsize, dptr);
    }
    pack_local_data_tm_.stop();
}

/* Merge data received from neighboring processor to local data. Append local
 * data by inserting new rows when necessary */
/* This code had to be modified to split existing row updates from new row
 * inserts for the sake of threading. Inserting a new row requires using a
 * critical section, which cannot be called within a branch. This modification
 * lets only master perform tasks related to inserting new rows. Current
 * implementation avoids the need for a critical section, but requires master to
 * loop over all rows and only insert what is needed.
 */
template <class T>
void PackedCommunicationBuffer::mergeRecvDataToMatrix(
    VariableSizeMatrix<T>& amat, const bool append)
{
    merge_or_insert_new_rows_tm_.start();

    if (append)
    {
        /* update existing rows and insert new rows. This part is written
         * this way to make it amenable to threading
         */
        std::vector<int> row_size(packed_num_rows_, 0);
        /* update existing entries*/
// Threading here leads to results slightly dependent on number of threads (jlf,
// 07/15/2016)
#pragma omp parallel for
        for (int i = 0; i < packed_num_rows_; i++)
        {
            int* rindex = (int*)amat.getTableValue(packed_lvars_ptr_[i]);
            const int k = packed_nnzrow_ptr_[i + 1] - packed_nnzrow_ptr_[i];
            row_size[i] = k;

            if (rindex != nullptr)
            {
                const int lrindex = *rindex;
                /* insert columns */
                amat.updateLocalRowAdd(k, lrindex,
                    &packed_pj_ptr_[packed_nnzrow_ptr_[i]],
                    &packed_pa_ptr_[packed_nnzrow_ptr_[i]]);
                row_size[i] = -1;
            }
        }

        /* now insert new rows */
        for (int i = 0; i < packed_num_rows_; i++)
        {
            if (row_size[i] != -1)
            {
                amat.insertNewRow(row_size[i], packed_lvars_ptr_[i],
                    &packed_pj_ptr_[packed_nnzrow_ptr_[i]],
                    &packed_pa_ptr_[packed_nnzrow_ptr_[i]], append);
            }
        }
    }
    else // append is false so update only existing rows
    {
        updateMatrixEntriesWithRecvBuf(amat);
    }

    merge_or_insert_new_rows_tm_.stop();
}

/* Update local data with data received from neighboring processor - Merge only
 * on existing data. No new data is merged */
template <class T>
void PackedCommunicationBuffer::updateMatrixEntriesWithRecvBuf(
    VariableSizeMatrix<T>& amat)
{
    merge_to_existing_rows_tm_.start();
#pragma omp parallel for
    for (int i = 0; i < packed_num_rows_; i++)
    {
        int* rindex = (int*)amat.getTableValue(packed_lvars_ptr_[i]);
        /* row on this proc. */
        if (rindex != nullptr)
        {
            const int lrindex = *rindex;
            /* insert columns */
            int k = packed_nnzrow_ptr_[i + 1] - packed_nnzrow_ptr_[i];
            amat.updateLocalRowAdd(k, lrindex,
                &packed_pj_ptr_[packed_nnzrow_ptr_[i]],
                &packed_pa_ptr_[packed_nnzrow_ptr_[i]]);
        }
    }
    merge_to_existing_rows_tm_.stop();
}

/* Copy rows from recv buffer and insert into matrix */
template <class T>
void PackedCommunicationBuffer::insertRowsFromRecvBuf(
    VariableSizeMatrix<T>& amat, const bool append)
{
    copy_and_insert_new_rows_tm_.start();
    for (int i = 0; i < packed_num_rows_; i++)
    {
        int* rindex     = (int*)amat.getTableValue(packed_lvars_ptr_[i]);
        const int start = packed_nnzrow_ptr_[i];
        const int ncols = packed_nnzrow_ptr_[i + 1] - start;
        /* row on this proc. */
        if (rindex != nullptr)
        {
            /* insert columns */
            amat.initializeLocalRow(
                ncols, *rindex, &packed_pj_ptr_[start], &packed_pa_ptr_[start]);
        }
        else
        {
            /* Insert new row */
            amat.insertNewRow(ncols, packed_lvars_ptr_[i],
                &packed_pj_ptr_[packed_nnzrow_ptr_[i]],
                &packed_pa_ptr_[packed_nnzrow_ptr_[i]], append);
        }
    }
    copy_and_insert_new_rows_tm_.stop();
}

template void PackedCommunicationBuffer::mergeRecvDataToMatrix(
    VariableSizeMatrix<SparseRow>& amat, const char* data, const bool append);
template void PackedCommunicationBuffer::mergeRecvDataToMatrix(
    VariableSizeMatrix<SparseRowAndTable>& amat, const char* data,
    const bool append);
template void PackedCommunicationBuffer::updateMatrixEntriesWithRecvBuf(
    VariableSizeMatrix<SparseRow>& amat, const char* data);
template void PackedCommunicationBuffer::updateMatrixEntriesWithRecvBuf(
    VariableSizeMatrix<SparseRowAndTable>& amat, const char* data);
template void PackedCommunicationBuffer::insertRowsFromRecvBuf(
    VariableSizeMatrix<SparseRow>& amat, const char* data, const bool append);
template void PackedCommunicationBuffer::insertRowsFromRecvBuf(
    VariableSizeMatrix<SparseRowAndTable>& amat, const char* data,
    const bool append);

template void PackedCommunicationBuffer::mergeRecvDataToMatrix(
    VariableSizeMatrix<SparseRow>& amat, const bool append);
template void PackedCommunicationBuffer::mergeRecvDataToMatrix(
    VariableSizeMatrix<SparseRowAndTable>& amat, const bool append);
template void PackedCommunicationBuffer::updateMatrixEntriesWithRecvBuf(
    VariableSizeMatrix<SparseRow>& amat);
template void PackedCommunicationBuffer::updateMatrixEntriesWithRecvBuf(
    VariableSizeMatrix<SparseRowAndTable>& amat);
template void PackedCommunicationBuffer::insertRowsFromRecvBuf(
    VariableSizeMatrix<SparseRow>& amat, const bool append);
template void PackedCommunicationBuffer::insertRowsFromRecvBuf(
    VariableSizeMatrix<SparseRowAndTable>& amat, const bool append);
