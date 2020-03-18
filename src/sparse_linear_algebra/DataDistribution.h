// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DATADISTRIBUTION_H_
#define MGMOL_DATADISTRIBUTION_H_

#include "DirectionalReduce.h"
#include "MPIdata.h"
#include "PEenv.h"
#include "VariableSizeMatrix.h"

#include <mpi.h>

#include <iostream>

class DataDistribution
{
    static Timer gathersizes_tm_;
    static Timer reducesizes_tm_;
    static Timer distribute_local_row_tm_;
    static Timer distribute_local_data_tm_;
    static Timer send_recv_tm_;
    static Timer send_recv_rows_tm_;
    static Timer update_row_tm_;
    static Timer augment_local_data_tm_;
    static Timer update_local_rows_tm_;
    static Timer initLocalRow_tm_;
    static Timer send_recv_ovlp_tm_;
    static Timer send_recv_rows_ovlp_tm_;
    static Timer send_recv_ovlp_wait_tm_;
    static Timer send_recv_rows_ovlp_wait_tm_;

    static int max_matsize_;
    static int max_nnz_;

    static short count_computeMaxDataSize_;
    static short maxcount_computeMaxDataSize_;

    std::string name_;
    const double spread_radius_; // Spreading radius for data distribution

    DirectionalReduce* dir_reduce_;

    const pb::PEenv& mypeenv_;

    int lsize_; // Initial size of local matrix
    int aug_size_; // Augmented size of the local matrix after data distribution
    MPI_Comm cart_comm_; // MPI cartesian communicator for data distribution
    int data_pos_[6]; // positions of datatypes within buffer - 2 entries per
                      // direction
    int loc_data_sz_[3]; // size of initial local data that is packed for data
                         // distribution

    /* pointers to communication recv buffer entries */
    /* position of number of rows of matrix in recv buffer */
    int* rbuf_nrows_ptr_;
    /* start position of double datatypes in buffer */
    int* rbuf_start_double_pos_ptr_;
    /* start position of matrix nonzero row count info */
    int* rbuf_nnzrow_ptr_;
    /* start position of local variables global indexes array info */
    int* rbuf_lvars_ptr_;
    /* start position of column index array info */
    int* rbuf_pj_ptr_;
    /* start position of matrix (double) coefficient data */
    double* rbuf_pa_ptr_;
    /* actual size of data in recv buffer (in bytes or sizeof char) */
    int rbuf_data_size_;

    //  template <class T>
    /* compute starting positions for packing local data */
    //  template <class T>
    void computePackedDataPositions(
        const VariableSizeMatrix<sparserow>& lmat, int* pos) const
    {
        const int lsize  = lmat.n();
        const int nnzmat = lmat.nnzmat();

        /* now compute start positions for packing local data */
        pos[0] = 0; /* position of integer variables */
        const int offset
            = 2 * (lsize + 1) * sizeof(int) + (nnzmat + 1) * sizeof(int);
        const int padding = (offset % sizeof(double));
        pos[1] = offset + padding; /* position of double variables - adjusted
                                      for alignment */

        return;
    }
    /* get maximum buffer size for data transfer */
    int getPackedBufferSize(const int maxsize, const int nzmax) const
    {
        /* compute buffer size for packing data - include alignment padding for
         * double, if any */
        const int offset
            = 2 * (maxsize + 1) * sizeof(int) + (nzmax + 1) * sizeof(int);
        const int padding = (offset % sizeof(double));
        const int bsiz    = offset + nzmax * sizeof(double) + padding;

        return bsiz;
    }

    /* Set data pointer positions on recv buffer */
    void setPointersToRecvData(const char* rbuf);

    /* Reset data pointer positions on recv buffer to NULL*/
    void resetPointersToRecvDataToNULL()
    {
        rbuf_nrows_ptr_            = nullptr;
        rbuf_start_double_pos_ptr_ = nullptr;
        rbuf_nnzrow_ptr_           = nullptr;
        rbuf_lvars_ptr_            = nullptr;
        rbuf_pj_ptr_               = nullptr;
        rbuf_pa_ptr_               = nullptr;

        rbuf_data_size_ = 0;
    }

    /* Perform data distribution of local data */
    template <class T>
    void distributeLocalData(const int nsteps, const int dir, const int disp,
        const int bsiz, const int* pos, VariableSizeMatrix<sparserow>& lmat,
        VariableSizeMatrix<T>& aug_mat, const bool append, const bool bcflag);
    template <class T>
    void distributeLocalDataWithCommOvlp(const int nsteps, const int dir,
        const int disp, const int bsiz, const int* pos,
        VariableSizeMatrix<sparserow>& lmat, VariableSizeMatrix<T>& aug_mat,
        const bool append, const bool bcflag);
    template <class T>
    void distributeLocalRows(const int nsteps, const int dir, const int disp,
        const int bsiz, const int* pos, VariableSizeMatrix<sparserow>& lmat,
        VariableSizeMatrix<T>& aug_mat, const bool append);
    template <class T>
    void distributeLocalRowsWithCommOvlp(const int nsteps, const int dir,
        const int disp, const int bsiz, const int* pos,
        VariableSizeMatrix<sparserow>& lmat, VariableSizeMatrix<T>& aug_mat,
        const bool append);
    template <class T>
    void mergeDataFromNeighborToLocalData(
        VariableSizeMatrix<T>& aug_mat, const char* rbuf, const bool append);
    template <class T>
    void updateExistingLocalDataEntriesWithRecvBuf(
        VariableSizeMatrix<T>& amat, const char* rbuf);
    template <class T>
    void copyRowsFromRecvBuf(
        VariableSizeMatrix<T>& amat, const char* rbuf, const bool append);
    /* Pack local data into buffer */
    void packLocalData(
        VariableSizeMatrix<sparserow>& lmat, const int* pos, char* buf);
    /* get recv buffer data size */
    int getRecvBufDataSize() { return rbuf_data_size_; }
    void reduceDataSizes(const short dir,
        const VariableSizeMatrix<sparserow>& lmat, int* maxsize, int* nzmax);

    void computeMaxDataSize(const short dir,
        const VariableSizeMatrix<sparserow>& lmat, int* maxsize, int* nzmax);

public:
    DataDistribution(const std::string& name, const double s_radius,
        const pb::PEenv& myPEenv, const double domain[]);
    DataDistribution(const std::string& name, const int max_steps[3],
        const pb::PEenv& myPEenv, const double domain[]);
    ~DataDistribution() { delete dir_reduce_; }

    static void enforceComputeMaxDataSize()
    {
        count_computeMaxDataSize_ = maxcount_computeMaxDataSize_ - 1;
    }

    template <class T>
    void augmentLocalData(VariableSizeMatrix<T>& vsmat, const bool append,
        const bool bcflag = false); // augment the local matrix

    template <class T>
    void updateLocalRows(VariableSizeMatrix<T>& vsmat,
        const bool append = false); // augment the local matrix

    static void printTimers(std::ostream& os); // print timers

    void printStats()
    {
        if (onpe0)
        {
            std::cout << "spread_radius = " << spread_radius_
                      << " initial size = " << lsize_
                      << " augmented size = " << aug_size_ << std::endl;
            dir_reduce_->printStats(std::cout);
        }
    }

    /* Communicate and fill matrix data.
     * Communicate data corresponding to gids and
     * receive and assemble data from neighboring procs.
     */
    template <class T>
    void consolidateMatrix(
        const std::vector<int>& gids, VariableSizeMatrix<T>& mat)
    {
        // reset matrix keeping only nonzero rows specified by gids
        mat.sparsify(gids);
        // gather/ distribute data from neighbors
        updateLocalRows(mat, true);
    }
};

#endif
