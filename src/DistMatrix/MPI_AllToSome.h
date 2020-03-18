// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: MPI_AllToSome.h 7 2011-03-08 18:50:31Z jeanluc $

template <typename ScalarType>
int MPI_AlltofirstN(ScalarType* sendbuf, const int blocksize,
    ScalarType* recvbuf, const int nfirst, MPI_Comm comm);

template <typename ScalarType>
int MPI_AlltofirstNv(ScalarType* sendbuf, int* sendcnts, int* sdispls,
    ScalarType* recvbuf, int* recvcnts, int* rdispls, const int nfirst,
    MPI_Comm comm);
