// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SymmetricMatrix.h"

template <>
void SymmetricMatrix<short>::mpiAllOr()
{
    short* recvbuf = new short[size_];
    MPI_Allreduce(&data_[0], recvbuf, size_, MPI_SHORT, MPI_BOR, comm_);
    memcpy(&data_[0], recvbuf, size_ * sizeof(short));

    delete[] recvbuf;
}

template <>
void SymmetricMatrix<char>::mpiAllOr()
{
    char* recvbuf = new char[size_];
    MPI_Allreduce(&data_[0], recvbuf, size_, MPI_CHAR, MPI_BOR, comm_);
    memcpy(&data_[0], recvbuf, size_ * sizeof(char));

    delete[] recvbuf;
}

template <>
void SymmetricMatrix<int>::mpiAllOr()
{
    int* recvbuf = new int[size_];
    MPI_Allreduce(&data_[0], recvbuf, size_, MPI_INT, MPI_BOR, comm_);
    memcpy(&data_[0], recvbuf, size_ * sizeof(int));

    delete[] recvbuf;
}
