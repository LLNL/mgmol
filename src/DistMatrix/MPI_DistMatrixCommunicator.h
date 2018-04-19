// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MPI_DistMatrixCOMMUNICATOR_H
#define MPI_DistMatrixCOMMUNICATOR_H

#if USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif


namespace dist_matrix{

class MPI_DistMatrixCommunicator
{
private:
    MPI_Comm comm_;

public:
    MPI_DistMatrixCommunicator()
    {
        comm_=-1;
    }

    ~MPI_DistMatrixCommunicator()
    {
        assert( comm_!=-1 );
        MPI_Comm_free(&comm_);
    }
    
    MPI_Comm& comm()
    {
        return comm_;
    }

};

}

#endif
