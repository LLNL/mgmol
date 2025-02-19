// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DistributedIonicData.h"
#include "Mesh.h"
#include "tools.h"

#include <string.h>

DistributedIonicData::DistributedIonicData(
    const std::vector<std::string>& local_names,
    const std::vector<double>& local_data)
    : ion_names_(local_names), data_(local_data)
{
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();

    // save cartesian communicator info
    cart_comm_     = myPEenv.cart_comm();
    const int disp = -1;
    for (short dir = 0; dir < 3; dir++)
        MPI_Cart_shift(cart_comm_, dir, disp, &source_[dir], &dest_[dir]);
}

int DistributedIonicData::pack(char* cbuff, double* dbuff)
{
    std::vector<double>::iterator itf = data_.begin();
    double* dptr                      = dbuff;

    int idx = 0;
    for (std::vector<std::string>::iterator it = ion_names_.begin();
         it != ion_names_.end(); ++it)
    {
        std::string s(*it);
        FixedLengthString t;
        strncpy(t.mystring, s.c_str(), IonData_MaxStrLength - 1);
        memcpy(&cbuff[idx], t.mystring, IonData_MaxStrLength);
        idx += IonData_MaxStrLength;

        for (short i = 0; i < 3; i++)
        {
            *(dptr++) = *(itf++);
        }
    }

    return ion_names_.size();
}

void DistributedIonicData::unpack(char*& cptr, double*& dptr, const short ndata)
{
    for (short i = 0; i < ndata; i++)
    {
        // get name
        std::string name(cptr, IonData_MaxStrLength);
        stripLeadingAndTrailingBlanks(name);
        ion_names_.push_back(name);
        cptr += IonData_MaxStrLength;

        // get force
        for (short j = 0; j < 3; j++)
            data_.push_back(*(dptr++));
    }
}

// augment Data
void DistributedIonicData::augmentData(const int nsteps, const int dir,
    const int disp, const int locSize, const int maxLocSize,
    DistributedIonicData& data2send)
{
    // setup buffers for data transfer
    int csize              = IonData_MaxStrLength * locSize;
    const int cbuff_size_r = IonData_MaxStrLength * maxLocSize;
    char* cbuff            = new char[cbuff_size_r];
    char* cbuff_r          = new char[cbuff_size_r];

    int dsize              = 3 * locSize;
    const int dbuff_size_r = 3 * maxLocSize;
    double* dbuff          = new double[dbuff_size_r];
    double* dbuff_r        = new double[dbuff_size_r];

    int* ibuff   = new int[1];
    int* ibuff_r = new int[1];

    // pack data into buffers
    ibuff[0] = data2send.pack(cbuff, dbuff);
    assert(ibuff[0] >= 0);
    assert(ibuff[0] < 10000);

    // transfer data in multiple directions
    MPI_Request request1[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Request request2[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Request request3[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };

    int src;
    int dst;
    if (disp == -1)
    {
        src = source_[dir];
        dst = dest_[dir];
    }
    else
    {
        src = dest_[dir];
        dst = source_[dir];
    }

    int step = 0;
    while (step < nsteps)
    {
        // chars data
        MPI_Irecv(
            cbuff_r, cbuff_size_r, MPI_CHAR, src, 0, cart_comm_, &request1[0]);
        MPI_Isend(cbuff, csize, MPI_CHAR, dst, 0, cart_comm_, &request1[1]);
        // doubles data
        MPI_Irecv(dbuff_r, dbuff_size_r, MPI_DOUBLE, src, 0, cart_comm_,
            &request2[0]);
        MPI_Isend(dbuff, dsize, MPI_DOUBLE, dst, 0, cart_comm_, &request2[1]);
        // ints data
        MPI_Irecv(ibuff_r, 1, MPI_INT, src, 0, cart_comm_, &request3[0]);
        MPI_Isend(ibuff, 1, MPI_INT, dst, 0, cart_comm_, &request3[1]);

        // wait to complete data transfer
        MPI_Waitall(2, request1, MPI_STATUS_IGNORE);
        MPI_Waitall(2, request2, MPI_STATUS_IGNORE);
        MPI_Waitall(2, request3, MPI_STATUS_IGNORE);

        // unpack: update ions_data for next call to this function
        char* cptr   = &cbuff_r[0];
        double* dptr = &dbuff_r[0];
        // first get data size
        const int rsize = *ibuff_r;
        assert(rsize >= 0);
        assert(rsize < 10000);

        // unpack received data and add it to data
        unpack(cptr, dptr, rsize);

        // done with this phase of data transfer

        // copy recv buffer to send buffer
        memcpy(cbuff, cbuff_r, rsize * IonData_MaxStrLength * sizeof(char));
        memcpy(dbuff, dbuff_r, rsize * 3 * sizeof(double));

        step++;
    }

    delete[] cbuff;
    delete[] cbuff_r;
    delete[] dbuff;
    delete[] dbuff_r;
    delete[] ibuff;
    delete[] ibuff_r;
}
