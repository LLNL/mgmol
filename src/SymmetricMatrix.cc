// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "SymmetricMatrix.h"
#include <string.h>

template <>
void SymmetricMatrix<short>::mpiAllOr()
{
#ifdef USE_MPI
    short* recvbuf=new short[size_];
    MPI_Allreduce(&data_[0], recvbuf, size_, 
                  MPI_SHORT, MPI_BOR, comm_);
    memcpy(&data_[0], recvbuf, size_*sizeof(short) );

    delete[] recvbuf;
#endif
}

template <>
void SymmetricMatrix<char>::mpiAllOr()
{
#ifdef USE_MPI
    char* recvbuf=new char[size_];
    MPI_Allreduce(&data_[0], recvbuf, size_, 
                  MPI_CHAR, MPI_BOR, comm_);
    memcpy(&data_[0], recvbuf, size_*sizeof(char) );

    delete[] recvbuf;
#endif
}

template <>
void SymmetricMatrix<int>::mpiAllOr()
{
#ifdef USE_MPI
    int* recvbuf=new int[size_];
    MPI_Allreduce(&data_[0], recvbuf, size_, 
                  MPI_INT, MPI_BOR, comm_);
    memcpy(&data_[0], recvbuf, size_*sizeof(int) );

    delete[] recvbuf;
#endif
}
