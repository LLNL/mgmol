// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef _PACKEDCOMMUNICATIONBUFFER_H_
#define _PACKEDCOMMUNICATIONBUFFER_H_

#include "VariableSizeMatrix.h"

#include <mpi.h>

class PackedCommunicationBuffer
{
    static Timer pack_local_data_tm_;
    static Timer merge_or_insert_new_rows_tm_;
    static Timer merge_to_existing_rows_tm_;
    static Timer copy_and_insert_new_rows_tm_;

    /* communication buffers */
    /* actual data buffer */
    static char* storage_;
    /* pointer to send buffer */
    static char* send_buffer_;
    /* pointer to recv buffer */
    static char* recv_buffer_;

    /* size of storage buffer */
    static int storage_size_;

    /* size of data in packed buffer (in bytes or sizeof char)*/
    int packed_data_size_;
    /* number of rows in packed buffer */
    int packed_num_rows_;
    /* start position of matrix nonzero row count info */
    int* packed_nnzrow_ptr_;
    /* start position of local variables global indexes array info */
    int* packed_lvars_ptr_;
    /* start position of column index array info */
    int* packed_pj_ptr_;
    /* start position of matrix (double) coefficient data */
    double* packed_pa_ptr_;

    /* Pack local data into buffer */
    void packLocalData(VariableSizeMatrix<sparserow>& lmat, const int* pos);

public:
    PackedCommunicationBuffer(const int bsize);
    ~PackedCommunicationBuffer()
    {
        // we don't delete storage_ here since it is static
        // and will be reused in other objects
    }

    static void deleteStorage()
    {
        if (storage_ != nullptr) delete[] storage_;
        storage_ = nullptr;
    }

    /* initialize send buffer with local matrix */
    void initialize(
        VariableSizeMatrix<sparserow>& lmat, const int* data_type_pos)
    {
        assert(storage_ != 0);
        packLocalData(lmat, data_type_pos);
    }

    char* sendBuffer() { return send_buffer_; }
    char* recvBuffer() { return recv_buffer_; }

    void swapSendRecvBuffers()
    {
        /* swap pointers for next communication step */
        char* tmp    = send_buffer_;
        send_buffer_ = recv_buffer_;
        recv_buffer_ = tmp;
    }

    /* get recv buffer data size */
    int getRecvBufDataSize() { return packed_data_size_; }

    /* setup pointers for data accumulation */
    void setupPackedDataPointers(const char* data);

    int getRecvBufDataSize(const char* data)
    {
        setupPackedDataPointers(data);
        return getRecvBufDataSize();
    }

    /* data assembly routines */
    /* Merge data received from neighboring processor to local data. Append
     * local data by inserting new rows when necessary */
    template <class T>
    void mergeRecvDataToMatrix(
        VariableSizeMatrix<T>& amat, const char* data, const bool append)
    {
        /* setup pointers for data accumulation into matrix amat */
        setupPackedDataPointers(data);
        mergeRecvDataToMatrix(amat, append);
    }
    /* Update local data with data received from neighboring processor - Merge
     * only on existing data. No new data is merged */
    template <class T>
    void updateMatrixEntriesWithRecvBuf(
        VariableSizeMatrix<T>& amat, const char* data)
    {
        /* setup pointers for data accumulation into matrix amat */
        setupPackedDataPointers(data);
        updateMatrixEntriesWithRecvBuf(amat);
    }
    /* Copy rows from recv buffer and insert into matrix */
    template <class T>
    void insertRowsFromRecvBuf(
        VariableSizeMatrix<T>& amat, const char* data, const bool append)
    {
        /* setup pointers for data accumulation into matrix amat */
        setupPackedDataPointers(data);
        insertRowsFromRecvBuf(amat, append);
    }

    template <class T>
    void mergeRecvDataToMatrix(
        VariableSizeMatrix<T>& aug_mat, const bool append);
    template <class T>
    void updateMatrixEntriesWithRecvBuf(VariableSizeMatrix<T>& amat);
    template <class T>
    void insertRowsFromRecvBuf(VariableSizeMatrix<T>& amat, const bool append);

    static void printTimers(std::ostream& os)
    {
        pack_local_data_tm_.print(os);
        // setup_data_ptr_tm_.print(os);
        merge_or_insert_new_rows_tm_.print(os);
        merge_to_existing_rows_tm_.print(os);
        copy_and_insert_new_rows_tm_.print(os);
    }
};

#endif
