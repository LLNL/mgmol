// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: MPI_AllToSome.h 7 2011-03-08 18:50:31Z jeanluc $
int MPI_AlltofirstN( double *sendbuf, const int blocksize,
                     double *recvbuf, const int nfirst,
                     MPI_Comm comm );
int MPI_AlltofirstN( int *sendbuf, const int blocksize,
                     int *recvbuf, const int nfirst,
                     MPI_Comm comm );
int MPI_AlltofirstN( unsigned short *sendbuf, const int blocksize,
                     unsigned short *recvbuf, const int nfirst,
                     MPI_Comm comm );

int MPI_AlltofirstNv( int *sendbuf, int *sendcnts, int *sdispls,
                      int *recvbuf, int *recvcnts, int *rdispls,
                      const int nfirst,
                      MPI_Comm comm );
int MPI_AlltofirstNv( unsigned short *sendbuf, int *sendcnts, int *sdispls,
                      unsigned short *recvbuf, int *recvcnts, int *rdispls,
                      const int nfirst,
                      MPI_Comm comm );
