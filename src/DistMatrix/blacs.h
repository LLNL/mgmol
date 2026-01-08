// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_BLACS_H
#define MGMOL_BLACS_H

#include <mpi.h>

extern "C"
{
    void igesd2d(int*, int*, int*, int*, int*, int*, int*);
    void sgesd2d(int*, int*, int*, double*, int*, int*, int*);
    void igerv2d(int*, int*, int*, int*, int*, int*, int*);
    void sgerv2d(int*, int*, int*, double*, int*, int*, int*);
    void sgsum2d(int*, char*, char*, int*, int*, double*, int*, int*, int*);
    void igamn2d(int*, char*, char*, int*, int*, int*, int*, int*, int*, int*,
        int*, int*);
    void blacs_pinfo(int*, int*);
    void blacs_get(int*, int*, int*);
    void blacs_barrier(int*, char*);
    void blacs_gridinfo(int*, int*, int*, int*, int*);
    void blacs_gridinit(int*, char*, int*, int*);
    void blacs_gridmap(int*, int*, int*, int*, int*);
    void blacs_abort(int*, int*);
    void blacs_gridexit(int*);
    void blacs_exit(int*);
    int blacs_pnum(int*, int*, int*);
}

#ifdef MGMOL_USE_SCALAPACK
extern "C"
{
#endif
    // C interface to the BLACS
    void Cdgesd2d(int, int, int, double*, int, int, int);
    void Cdgerv2d(int, int, int, double*, int, int, int);
    void Cdgsum2d(int, char*, char*, int, int, double*, int, int, int);

    void Cblacs_pinfo(int*, int*);
    void Cblacs_get(int, int, int*);
    void Cblacs_barrier(int, char*);
    void Cblacs_gridinfo(int, int*, int*, int*, int*);
    void Cblacs_gridinit(int*, char[], int, int);
    void Cblacs_gridmap(int*, int*, int, int, int);
    void Cblacs_abort(int, int);
    void Cblacs_gridexit(int);
    void Cblacs_exit(int);
    int Cblacs_pnum(int, int, int);
    int Csys2blacs_handle(MPI_Comm);
#ifdef MGMOL_USE_SCALAPACK
}
#endif

#endif
