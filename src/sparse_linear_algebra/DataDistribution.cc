// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DataDistribution.h"
#include "../tools.h"
#include "PackedCommunicationBuffer.h"
#include "VariableSizeMatrix.h"

#include "../Control.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <mpi.h>

Timer DataDistribution::gathersizes_tm_("DataDistribution::gatherDataSizes");
Timer DataDistribution::reducesizes_tm_("DataDistribution::reduceDataSizes");
Timer DataDistribution::distribute_local_row_tm_(
    "DataDistribution::distributeLocalRows");
Timer DataDistribution::distribute_local_data_tm_(
    "DataDistribution::distributeLocalData");
Timer DataDistribution::send_recv_tm_("DataDistribution::MPISendRecv");
Timer DataDistribution::send_recv_rows_tm_(
    "DataDistribution::MPISendRecv_rows");
// Timer
// DataDistribution::update_row_tm_("DataDistribution::consolidateLocalRows");
Timer DataDistribution::augment_local_data_tm_(
    "DataDistribution::augmentLocalData");
Timer DataDistribution::update_local_rows_tm_(
    "DataDistribution::updateLocalRows");
Timer DataDistribution::initLocalRow_tm_("DataDistribution::initLocalRows");
Timer DataDistribution::send_recv_ovlp_tm_(
    "DataDistribution::MPISendRecv_ovlp");
Timer DataDistribution::send_recv_ovlp_wait_tm_(
    "DataDistribution::MPISendRecv_ovlp_wait");
Timer DataDistribution::send_recv_rows_ovlp_tm_(
    "DataDistribution::MPISendRecv_Rows_ovlp");
Timer DataDistribution::send_recv_rows_ovlp_wait_tm_(
    "DataDistribution::MPISendRecv_Rows_ovlp_wait");

int DataDistribution::max_matsize_                   = -1;
int DataDistribution::max_nnz_                       = -1;
short DataDistribution::count_computeMaxDataSize_    = 0;
short DataDistribution::maxcount_computeMaxDataSize_ = 500;

/* Constructor */
DataDistribution::DataDistribution(const std::string& name,
    const int max_steps[3], const pb::PEenv& myPEenv, const double[])
    : name_(name), spread_radius_(0.0), mypeenv_(myPEenv)
{
    /* get cartesian communicator */
    cart_comm_ = myPEenv.cart_comm();

    dir_reduce_ = new DirectionalReduce(myPEenv.cart_comm(), max_steps);

    rbuf_nrows_ptr_            = nullptr;
    rbuf_start_double_pos_ptr_ = nullptr;
    rbuf_nnzrow_ptr_           = nullptr;
    rbuf_lvars_ptr_            = nullptr;
    rbuf_pj_ptr_               = nullptr;
    rbuf_pa_ptr_               = nullptr;
    rbuf_data_size_            = 0;

    lsize_ = -1;
}

/* Constructor */
DataDistribution::DataDistribution(const std::string& name,
    const double s_radius, const pb::PEenv& myPEenv, const double domain[])
    : name_(name), spread_radius_(s_radius), mypeenv_(myPEenv)
{
    assert(domain[0] > 0.);
    assert(domain[1] > 0.);
    assert(domain[2] > 0.);
    assert(s_radius > 0.);

    /* get cartesian communicator */
    cart_comm_ = myPEenv.cart_comm();

    dir_reduce_ = new DirectionalReduce(myPEenv.cart_comm(), s_radius, domain);

    assert(spread_radius_ > 0.0);

    rbuf_nrows_ptr_            = nullptr;
    rbuf_start_double_pos_ptr_ = nullptr;
    rbuf_nnzrow_ptr_           = nullptr;
    rbuf_lvars_ptr_            = nullptr;
    rbuf_pj_ptr_               = nullptr;
    rbuf_pa_ptr_               = nullptr;
    rbuf_data_size_            = 0;

    lsize_ = -1;
}

/* Merge data received from neighboring processor to local data. Append local
 * data by inserting new rows when necessary */
template <class T>
void DataDistribution::mergeDataFromNeighborToLocalData(
    VariableSizeMatrix<T>& amat, const char* rbuf, const bool append)
{
    //   const int n = *rbuf_nrows_ptr_;
    //   int* const lvars = rbuf_lvars_ptr_;
    //   int* const nnzrow = rbuf_nnzrow_ptr_;
    //   int* const pj = rbuf_pj_ptr_;
    //   double* const pa = rbuf_pa_ptr_;

    /* set pointers to data in rbuf */
    setPointersToRecvData(rbuf);

    for (int i = 0; i < *rbuf_nrows_ptr_; i++)
    {
        int* rindex = (int*)amat.getTableValue(rbuf_lvars_ptr_[i]);
        const int k = rbuf_nnzrow_ptr_[i + 1] - rbuf_nnzrow_ptr_[i];
        /* row on this proc. */
        if (rindex != nullptr)
        {
            const int lrindex = *rindex;
            /* insert columns */
            amat.updateLocalRowAdd(k, lrindex,
                &rbuf_pj_ptr_[rbuf_nnzrow_ptr_[i]],
                &rbuf_pa_ptr_[rbuf_nnzrow_ptr_[i]]);
        }
        else
        {
            /* Insert new row */
            amat.insertNewRow(k, rbuf_lvars_ptr_[i],
                &rbuf_pj_ptr_[rbuf_nnzrow_ptr_[i]],
                &rbuf_pa_ptr_[rbuf_nnzrow_ptr_[i]], append);
        }
    }
}

/* Update local data with data received from neighboring processor - Merge only
 * on existing data. No new data is merged */
template <class T>
void DataDistribution::updateExistingLocalDataEntriesWithRecvBuf(
    VariableSizeMatrix<T>& amat, const char* rbuf)
{
    //   const int n = *rbuf_nrows_ptr_;
    //   int* const lvars = rbuf_lvars_ptr_;
    //   int* const nnzrow = rbuf_nnzrow_ptr_;
    //   int* const pj = rbuf_pj_ptr_;
    //   double* const pa = rbuf_pa_ptr_;

    /* set pointers to data in rbuf */
    setPointersToRecvData(rbuf);

    for (int i = 0; i < *rbuf_nrows_ptr_; i++)
    {
        int* rindex = (int*)amat.getTableValue(rbuf_lvars_ptr_[i]);
        /* row on this proc. */
        if (rindex != nullptr)
        {
            const int lrindex = *rindex;
            /* insert columns */
            int k = rbuf_nnzrow_ptr_[i + 1] - rbuf_nnzrow_ptr_[i];
            amat.updateLocalRowAdd(k, lrindex,
                &rbuf_pj_ptr_[rbuf_nnzrow_ptr_[i]],
                &rbuf_pa_ptr_[rbuf_nnzrow_ptr_[i]]);
        }
    }
}

template <class T>
void DataDistribution::copyRowsFromRecvBuf(
    VariableSizeMatrix<T>& amat, const char* const rbuf, const bool append)
{
    //   const int n = *rbuf_nrows_ptr_;
    //   int* const lvars = rbuf_lvars_ptr_;
    //   int* const nnzrow = rbuf_nnzrow_ptr_;
    //   int* const pj = rbuf_pj_ptr_;
    //   double* const pa = rbuf_pa_ptr_;

    /* set pointers to data in rbuf */
    setPointersToRecvData(rbuf);

    initLocalRow_tm_.start();
    for (int i = 0; i < *rbuf_nrows_ptr_; i++)
    {
        int* rindex     = (int*)amat.getTableValue(rbuf_lvars_ptr_[i]);
        const int start = rbuf_nnzrow_ptr_[i];
        const int ncols = rbuf_nnzrow_ptr_[i + 1] - start;
        /* row on this proc. */
        if (rindex != nullptr)
        {
            /* insert columns */
            amat.initializeLocalRow(
                ncols, *rindex, &rbuf_pj_ptr_[start], &rbuf_pa_ptr_[start]);
        }
        else
        {
            /* Insert new row */
            amat.insertNewRow(ncols, rbuf_lvars_ptr_[i],
                &rbuf_pj_ptr_[rbuf_nnzrow_ptr_[i]],
                &rbuf_pa_ptr_[rbuf_nnzrow_ptr_[i]], append);
        }
    }
    initLocalRow_tm_.stop();
}

/* Pack local data into buffer */
// template <class T>
void DataDistribution::packLocalData(
    VariableSizeMatrix<sparserow>& lmat, const int* pos, char* buf)
{

    /* Get some data info */
    const int lsize = lmat.n();

    /* Manually pack buffer */
    /* check for zero data - base case */
    if (lsize == 0)
    {
        int* const iptr = (int*)&buf[pos[0]];
        *(iptr)         = 0;
    }
    else
    {
        /* pack matrix size - n */
        int* iptr = (int*)&buf[pos[0]];
        *(iptr++) = lsize;
        /* pack start position of doubles (matrix entries) */
        *(iptr++) = (int)pos[1];
        /* pack nnzrow */
        *(iptr++) = 0;
        int k     = 0;
        for (int i = 0; i < lsize; i++)
        {
            k += lmat.nnzrow(i);
            *(iptr++) = k;
        }
        // pack matrix data
        double* dptr = (double*)&buf[pos[1]];
        lmat.copyDataToArray(iptr, iptr + lsize, dptr);
    }
}

/* Perform data distribution of local data */
/* Data is packed in fixed-size CSR format */
template <class T>
void DataDistribution::distributeLocalDataWithCommOvlp(const int nsteps,
    const int dir, const int disp, const int bsiz, const int* pos,
    VariableSizeMatrix<sparserow>& lmat, VariableSizeMatrix<T>& amat,
    const bool append, const bool bcflag)
{
    // use barriers to ensure processors begin at the same time before timing
    // this routine
    //  MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    //  mmpi.barrier();
    distribute_local_data_tm_.start();

    PackedCommunicationBuffer packed_buffer(bsiz);

    /* initialize siz to size of local data */
    int siz = loc_data_sz_[dir];
    // if(siz > bsiz)cout<<"siz = "<<siz<<" and bsiz = "<<bsiz<<endl;
    /* check that buffer size is large enough */
    assert(siz <= bsiz);

    packed_buffer.initialize(lmat, pos);

    /* Begin data distribution */
    MPI_Request request[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };

    /* Get source and destination ID to send and recv data */
    int source, dest;
    MPI_Cart_shift(cart_comm_, dir, disp, &source, &dest);

    if (nsteps > 0)
    {
#ifndef NDEBUG // check if receive buffers are large enough
        // #if 1
        int remote_size;
        MPI_Irecv(&remote_size, 1, MPI_INT, source, 0, cart_comm_, &request[0]);
        MPI_Isend(&siz, 1, MPI_INT, dest, 0, cart_comm_, &request[1]);
        /* wait to complete communication */
        MPI_Waitall(2, request, MPI_STATUSES_IGNORE);
        if (remote_size > bsiz)
        {
            std::cout << "ERROR: " << name_ << ", dir=" << dir
                      << ", remote_size=" << remote_size << ", bsiz=" << bsiz
                      << std::endl;
            MPI_Abort(cart_comm_, EXIT_FAILURE);
        }
        // string stamp="DataDistribution ("+name_+"), buffer size checked...";
        // printWithTimeStamp(stamp,cout);
        MPI_Barrier(cart_comm_);
#endif
        // step 0: send local data
        /* Send and receive data */
        send_recv_tm_.start();
        int mpircv = MPI_Irecv(packed_buffer.recvBuffer(), bsiz, MPI_CHAR,
            source, 0, cart_comm_, &request[0]);
        if (mpircv != MPI_SUCCESS)
        {
            std::cout << "ERROR in MPI_Irecv, code=" << mpircv << std::endl;
            MPI_Abort(cart_comm_, EXIT_FAILURE);
        }
        int mpisnd = MPI_Isend(packed_buffer.sendBuffer(), siz, MPI_CHAR, dest,
            0, cart_comm_, &request[1]);
        if (mpisnd != MPI_SUCCESS)
        {
            std::cout << "ERROR in MPI_Isend, code=" << mpisnd << std::endl;
            MPI_Abort(cart_comm_, EXIT_FAILURE);
        }
        /* wait to complete communication */
        MPI_Waitall(2, request, MPI_STATUSES_IGNORE);
        send_recv_tm_.stop();

        // complete remaining data transfer steps -- overlap communication and
        // work
        for (int step = 1; step < nsteps; step++)
        {
            /* Prepare to send next data */
            /* set pointers to recv'd data to be accumulated later */
            packed_buffer.setupPackedDataPointers(packed_buffer.recvBuffer());
            /* get size of data in recv buffer */
            siz = packed_buffer.getRecvBufDataSize();
            /* check that buffer size is large enough */
            assert(siz <= bsiz);
            /* swap pointers for next communication step */
            packed_buffer.swapSendRecvBuffers();

            // post a send and recv
            send_recv_ovlp_tm_.start();
            MPI_Irecv(packed_buffer.recvBuffer(), bsiz, MPI_CHAR, source, 0,
                cart_comm_, &request[0]);
            MPI_Isend(packed_buffer.sendBuffer(), siz, MPI_CHAR, dest, 0,
                cart_comm_, &request[1]);
            send_recv_ovlp_tm_.stop();

            // merge previously recv'd data to local data -- now stored in
            // send_buffer uses previously set internal pointers to data to
            // accomplish this See PackedCommunicationBuffer class for details
            packed_buffer.mergeRecvDataToMatrix(amat, append);

            // wait to complete send/recv
            send_recv_ovlp_wait_tm_.start();
            MPI_Waitall(2, request, MPI_STATUSES_IGNORE);
            send_recv_ovlp_wait_tm_.stop();
        }
        // merge final recv'd data
        packed_buffer.setupPackedDataPointers(packed_buffer.recvBuffer());
        if (bcflag == true) // this is the last step - apply boundary condition
        {
            packed_buffer.updateMatrixEntriesWithRecvBuf(amat);
        }
        else
        {
            packed_buffer.mergeRecvDataToMatrix(amat, append);
        }
    }
    distribute_local_data_tm_.stop();
}

/* Perform data distribution of local data */
/* Data is packed in fixed-size CSR format */
template <class T>
void DataDistribution::distributeLocalData(const int nsteps, const int dir,
    const int disp, const int bsiz, const int* pos,
    VariableSizeMatrix<sparserow>& lmat, VariableSizeMatrix<T>& amat,
    const bool append, const bool bcflag)
{
    distribute_local_data_tm_.start();

    /* initialize siz to size of local data */
    int siz = loc_data_sz_[dir];
    /* check that buffer size is large enough */
    assert(siz <= bsiz);

    /* buffer for packing data for communication */
    PackedCommunicationBuffer packed_buffer(bsiz);
    /* pack local data */
    packed_buffer.initialize(lmat, pos);
    /* Begin data distribution */
    MPI_Request request[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };

    /* Get source and destination ID to send and recv data */
    int source, dest;
    MPI_Cart_shift(cart_comm_, dir, disp, &source, &dest);

    int step = 0;
    while (step < nsteps)
    {
        /* Send and receive data */
        send_recv_tm_.start();
        MPI_Irecv(packed_buffer.recvBuffer(), bsiz, MPI_CHAR, source, 0,
            cart_comm_, &request[0]);
        MPI_Isend(packed_buffer.sendBuffer(), siz, MPI_CHAR, dest, 0,
            cart_comm_, &request[1]);
        MPI_Waitall(2, request, MPI_STATUSES_IGNORE);
        send_recv_tm_.stop();

        if ((step == nsteps - 1)
            && (bcflag
                == true)) // this is the last step - apply boundary condition
        {
            packed_buffer.updateMatrixEntriesWithRecvBuf(
                amat, packed_buffer.recvBuffer());
        }
        else
        {
            packed_buffer.mergeRecvDataToMatrix(
                amat, packed_buffer.recvBuffer(), append);
        }

        /* Prepare to send next data. */
        siz = packed_buffer.getRecvBufDataSize();
        /* check that buffer size is large enough */
        assert(siz <= bsiz);
        /* swap pointers for next communication step */
        packed_buffer.swapSendRecvBuffers();

        step++;
    }
    distribute_local_data_tm_.stop();
}

/* Augment the local matrix - data distribution */
template <class T>
void DataDistribution::augmentLocalData(
    VariableSizeMatrix<T>& vsmat, const bool append, const bool bcflag)
{
    assert(dir_reduce_->lstep(0) >= 0);
    assert(dir_reduce_->lstep(1) >= 0);
    assert(dir_reduce_->lstep(2) >= 0);

    augment_local_data_tm_.start();

    lsize_ = vsmat.n();

    int maxsize, nzmax;
    // short has_datasize_converged=0;
    for (short dir = 0; dir < 3; dir++)
    {
        if (dir_reduce_->lstep(dir) > 0)
        {
            /* Spread and receive data in the x/y/z-direction.
             * This corresponds to spreading data
             * across the rows/cols/depth of the cartesian grid
             */
            /**************** PHASE I *******************/

            /* Get initial local data (or initial data to be sent) */
            VariableSizeMatrix<sparserow> mat(vsmat, false);
            computeMaxDataSize(dir, mat, &maxsize, &nzmax);

            /* setup some variables */
            int pos = 2 * dir;
            computePackedDataPositions(mat, &data_pos_[pos]);
            loc_data_sz_[dir]
                = (&data_pos_[pos])[1] + mat.nnzmat() * sizeof(double);

            const int buffer_size = getPackedBufferSize(maxsize, nzmax);

            /* determine whether or not to apply boundary condition */
            bool trim = false;
            if (dir_reduce_->lstep(dir) == dir_reduce_->rstep(dir))
                trim = bcflag;
            /* send data to the left and recv from right */
            short disp = -1;
            distributeLocalDataWithCommOvlp(dir_reduce_->lstep(dir), dir, disp,
                buffer_size, &data_pos_[pos], mat, vsmat, append, trim);
            //        distributeLocalData(lstep_[dir], dir, disp, buffer_size,
            //        &data_pos_[pos], mat, vsmat, append, trim);
            /**************** PHASE II *******************/

            /* send data to the right and recv from left */
            disp = 1;
            distributeLocalDataWithCommOvlp(dir_reduce_->rstep(dir), dir, disp,
                buffer_size, &data_pos_[pos], mat, vsmat, append, trim);
            //        distributeLocalData(rstep_[dir], dir, disp, buffer_size,
            //        &data_pos_[pos], mat, vsmat, append, trim);
        }
    }

    aug_size_ = vsmat.n();
    augment_local_data_tm_.stop();

    return;
}

/* Perform data distribution of local data */
/* Data is packed in fixed-size CSR format */
template <class T>
void DataDistribution::distributeLocalRowsWithCommOvlp(const int nsteps,
    const int dir, const int disp, const int bsiz, const int* pos,
    VariableSizeMatrix<sparserow>& lmat, VariableSizeMatrix<T>& amat,
    const bool append)
{
    // use barriers to ensure processors begin at the same time before timing
    // this routine
    //  MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    //  mmpi.barrier();
    distribute_local_row_tm_.start();

    /* initialize siz to size of local data */
    int siz = loc_data_sz_[dir];
    /* check that buffer size is large enough */
    assert(siz <= bsiz);

    /* buffer for packing data for communication */
    PackedCommunicationBuffer packed_buffer(bsiz);
    /* pack local data */
    packed_buffer.initialize(lmat, pos);
    /* Begin data distribution */
    MPI_Request request[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };

    /* Get source and destination ID to send and recv data */
    int source, dest;
    MPI_Cart_shift(cart_comm_, dir, disp, &source, &dest);

    if (nsteps > 0)
    {
        // step 0: send local data
        /* Send and receive data */
        send_recv_rows_tm_.start();
        MPI_Irecv(packed_buffer.recvBuffer(), bsiz, MPI_CHAR, source, 0,
            cart_comm_, &request[0]);
        MPI_Isend(packed_buffer.sendBuffer(), siz, MPI_CHAR, dest, 0,
            cart_comm_, &request[1]);
        /* wait to complete communication */
        MPI_Waitall(2, request, MPI_STATUSES_IGNORE);
        send_recv_rows_tm_.stop();

        // complete remaining data transfer steps -- overlap communication and
        // work
        for (int step = 1; step < nsteps; step++)
        {
            /* Prepare to send next data */
            /* set pointers to recv'd data to be accumulated later */
            packed_buffer.setupPackedDataPointers(packed_buffer.recvBuffer());
            /* get size of data in recv buffer */
            siz = packed_buffer.getRecvBufDataSize();
            /* check that buffer size is large enough */
            assert(siz <= bsiz);
            /* swap pointers for next communication step */
            packed_buffer.swapSendRecvBuffers();

            // post a send and recv
            send_recv_rows_ovlp_tm_.start();
            MPI_Irecv(packed_buffer.recvBuffer(), bsiz, MPI_CHAR, source, 0,
                cart_comm_, &request[0]);
            MPI_Isend(packed_buffer.sendBuffer(), siz, MPI_CHAR, dest, 0,
                cart_comm_, &request[1]);
            send_recv_rows_ovlp_tm_.stop();

            // copy previously recv'd data to local data -- now stored in
            // send_buffer uses previously set internal pointers to data to
            // accomplish this See PackedCommunicationBuffer class for details
            packed_buffer.insertRowsFromRecvBuf(amat, append);

            // wait to complete send/recv
            send_recv_rows_ovlp_wait_tm_.start();
            MPI_Waitall(2, request, MPI_STATUSES_IGNORE);
            send_recv_rows_ovlp_wait_tm_.stop();
        }
        // copy final recv'd data
        packed_buffer.setupPackedDataPointers(packed_buffer.recvBuffer());
        packed_buffer.insertRowsFromRecvBuf(amat, append);
    }
    distribute_local_row_tm_.stop();
}

/* Perform data distribution of local data */
/* Data is packed in fixed-size CSR format */
template <class T>
void DataDistribution::distributeLocalRows(const int nsteps, const int dir,
    const int disp, const int bsiz, const int* pos,
    VariableSizeMatrix<sparserow>& lmat, VariableSizeMatrix<T>& amat,
    const bool append)
{
    distribute_local_row_tm_.start();

    /* initialize siz to size of local data */
    int siz = loc_data_sz_[dir];
    /* check that buffer size is large enough */
    assert(siz <= bsiz);

    /* buffer for packing data for communication */
    PackedCommunicationBuffer packed_buffer(bsiz);
    /* pack local data */
    packed_buffer.initialize(lmat, pos);
    /* Begin data distribution */
    MPI_Request request[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };

    /* Get source and destination ID to send and recv data */
    int source, dest;
    MPI_Cart_shift(cart_comm_, dir, disp, &source, &dest);

    int step = 0;
    while (step < nsteps)
    {
        /* Send and receive data */
        send_recv_rows_tm_.start();
        MPI_Irecv(packed_buffer.recvBuffer(), bsiz, MPI_CHAR, source, 0,
            cart_comm_, &request[0]);
        MPI_Isend(packed_buffer.sendBuffer(), siz, MPI_CHAR, dest, 0,
            cart_comm_, &request[1]);
        MPI_Waitall(2, request, MPI_STATUSES_IGNORE);
        send_recv_rows_tm_.stop();

        /* copy data from recv buffer into local matrix */
        packed_buffer.insertRowsFromRecvBuf(
            amat, packed_buffer.recvBuffer(), append);

        /* Prepare to send next data. */
        siz = packed_buffer.getRecvBufDataSize();
        /* check that buffer size is large enough */
        assert(siz <= bsiz);
        /* swap pointers for next communication step */
        packed_buffer.swapSendRecvBuffers();

        step++;
    }
    distribute_local_row_tm_.stop();
}

/* Augment the local matrix - data distribution */
template <class T>
void DataDistribution::updateLocalRows(
    VariableSizeMatrix<T>& vsmat, const bool append)
{
    /* Local variables */
    short disp;
    lsize_ = vsmat.n();

    update_local_rows_tm_.start();

    int maxsize, nzmax;
    for (short dir = 0; dir < 3; dir++)
    {
        if (dir_reduce_->lstep(dir) > 0)
        {
            /* Spread and receive data in the x/y/z-direction.
             * This corresponds to spreading data
             * across the rows/cols/depth of the cartesian grid
             */
            /**************** PHASE I *******************/

            /* Get initial local data (or initial data to be sent) */
            VariableSizeMatrix<sparserow> mat(vsmat, false);
            computeMaxDataSize(dir, mat, &maxsize, &nzmax);

            /* setup some variables */
            short pos = 2 * dir;
            computePackedDataPositions(mat, &data_pos_[pos]);
            loc_data_sz_[dir]
                = (&data_pos_[pos])[1] + mat.nnzmat() * sizeof(double);
            const int buffer_size = getPackedBufferSize(maxsize, nzmax);

            /* send data to the left and recv from right */
            disp = -1;
            distributeLocalRowsWithCommOvlp(dir_reduce_->lstep(dir), dir, disp,
                buffer_size, &data_pos_[pos], mat, vsmat, append);
            //        distributeLocalRows(lstep_[dir], dir, disp, buffer_size,
            //        &data_pos_[pos], mat, vsmat, append);
            /**************** PHASE II *******************/

            /* send data to the right and recv from left */
            disp = 1;
            distributeLocalRowsWithCommOvlp(dir_reduce_->rstep(dir), dir, disp,
                buffer_size, &data_pos_[pos], mat, vsmat, append);
            //        distributeLocalRows(rstep_[dir], dir, disp, buffer_size,
            //        &data_pos_[pos], mat, vsmat, append);
        }
    }

    aug_size_ = vsmat.n();
    update_local_rows_tm_.stop();

    return;
}

void DataDistribution::computeMaxDataSize(const short dir,
    const VariableSizeMatrix<sparserow>& lmat, int* maxsize, int* nzmax)
{
    // limit number of calls to gatherDataSizes()
    // somewhat arbitrary and may need to be tuned
    // assumes buffer size needed doesn't increase too much after that

    if (count_computeMaxDataSize_ < maxcount_computeMaxDataSize_)
    {
        //      gatherDataSizes(dir,lmat,maxsize,nzmax);
        int data[2] = { lmat.n(), lmat.nnzmat() };
        dir_reduce_->computeDirMax(dir, data);
        *maxsize = data[0];
        *nzmax   = data[1];

        max_matsize_ = (*maxsize > max_matsize_) ? *maxsize : max_matsize_;
        max_nnz_     = (*nzmax > max_nnz_) ? *nzmax : max_nnz_;
    }

    if (count_computeMaxDataSize_ == maxcount_computeMaxDataSize_)
    {
        max_matsize_ *= 4.0; // 1.5; Changed by Ian to facilitate long O(N)
                             // runs.
        max_nnz_ *= 4.0; // 1.5; Changed by Ian to facilitate long O(N) runs.

        Control& ct = *(Control::instance());
        if (onpe0 && ct.verbose > 0)
            std::cout << "DataDistribution: Use maxsize=" << max_matsize_
                      << ", nzmax=" << max_nnz_ << std::endl;
    }

    *maxsize = max_matsize_;
    *nzmax   = max_nnz_;

    count_computeMaxDataSize_++;
}

// template <class T>
void DataDistribution::reduceDataSizes(const short dir,
    const VariableSizeMatrix<sparserow>& lmat, int* maxsize, int* nzmax)
{
    // use barriers to ensure processors begin at the same time before timing
    // this routine
    //  MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    //  mmpi.barrier();
    reducesizes_tm_.start();

    /* initialize send and recv buffer */
    int buf[2] = { lmat.n(), lmat.nnzmat() };

    if (dir == 0) mypeenv_.maxXdir(buf, 2);
    if (dir == 1) mypeenv_.maxYdir(buf, 2);
    if (dir == 2) mypeenv_.maxZdir(buf, 2);

    /* now assign values */
    *maxsize = buf[0];
    *nzmax   = buf[1];

    reducesizes_tm_.stop();

    return;
}
void DataDistribution::setPointersToRecvData(const char* rbuf)
{
    assert(rbuf != nullptr);

    rbuf_nrows_ptr_ = (int*)rbuf;
    /* check if matrix is nonzero */
    if (*rbuf_nrows_ptr_ != 0)
    {
        rbuf_start_double_pos_ptr_ = rbuf_nrows_ptr_ + 1;
        rbuf_nnzrow_ptr_           = rbuf_start_double_pos_ptr_ + 1;
        rbuf_lvars_ptr_            = rbuf_nnzrow_ptr_ + (*rbuf_nrows_ptr_ + 1);
        rbuf_pj_ptr_               = rbuf_lvars_ptr_ + (*rbuf_nrows_ptr_);
        rbuf_pa_ptr_ = (double*)(rbuf + *rbuf_start_double_pos_ptr_);

        /* compute size of data in buffer */
        rbuf_data_size_ = *rbuf_start_double_pos_ptr_
                          + rbuf_nnzrow_ptr_[*rbuf_nrows_ptr_] * sizeof(double);
    }
    else
    {
        /* The size is equal to the size of an int if nrows == 0. */
        rbuf_data_size_ = sizeof(int);
    }
}

void DataDistribution::printTimers(std::ostream& os)
{
    gathersizes_tm_.print(os);
    reducesizes_tm_.print(os);
    distribute_local_row_tm_.print(os);
    distribute_local_data_tm_.print(os);
    send_recv_tm_.print(os);
    send_recv_ovlp_tm_.print(os);
    send_recv_ovlp_wait_tm_.print(os);
    send_recv_rows_tm_.print(os);
    send_recv_rows_ovlp_tm_.print(os);
    send_recv_rows_ovlp_wait_tm_.print(os);
    //   update_row_tm_.print(os);
    augment_local_data_tm_.print(os);
    update_local_rows_tm_.print(os);
    initLocalRow_tm_.print(os);
}

template void DataDistribution::augmentLocalData(
    VariableSizeMatrix<SparseRow>& vsmat, const bool append, const bool bcflag);
template void DataDistribution::augmentLocalData(
    VariableSizeMatrix<SparseRowAndTable>& vsmat, const bool append,
    const bool bcflag);
template void DataDistribution::updateLocalRows(
    VariableSizeMatrix<SparseRow>& vsmat, const bool append);
template void DataDistribution::updateLocalRows(
    VariableSizeMatrix<SparseRowAndTable>& vsmat, const bool append);
template void DataDistribution::consolidateMatrix(
    const std::vector<int>& gids, VariableSizeMatrix<SparseRow>& mat);
template void DataDistribution::consolidateMatrix(
    const std::vector<int>& gids, VariableSizeMatrix<SparseRowAndTable>& mat);
