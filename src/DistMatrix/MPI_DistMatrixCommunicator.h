// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MPI_DistMatrixCOMMUNICATOR_H
#define MPI_DistMatrixCOMMUNICATOR_H

#include <mpi.h>

namespace dist_matrix
{

class MPI_DistMatrixCommunicator
{
private:
    MPI_Comm comm_;

public:
    MPI_DistMatrixCommunicator() { comm_ = MPI_COMM_NULL; }

    ~MPI_DistMatrixCommunicator()
    {
        assert(comm_ != MPI_COMM_NULL);
        MPI_Comm_free(&comm_);
    }

    MPI_Comm& comm() { return comm_; }
};
}

#endif
