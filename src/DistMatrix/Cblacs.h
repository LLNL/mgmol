// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Cblacs.h 7 2011-03-08 18:50:31Z jeanluc $
#ifndef CBLACS_H
#define CBLACS_H

#ifdef Add_
#define sgerv2d sgerv2d_
#define dgsum2d dgsum2d_
#endif

#ifdef CRAY_T3E

#include <fortran.h>
#define CCHAR _fcd

#else

#define CCHAR char*

#endif

#ifndef UpCase

#define BLACS_BARRIER blacs_barrier
#define BLACS_PINFO blacs_pinfo
#define BLACS_GET blacs_get
#define BLACS_GRIDINIT blacs_gridinit
#define BLACS_GRIDMAP blacs_gridmap
#define BLACS_ABORT blacs_abort
#define BLACS_GRIDEXIT blacs_gridexit
#define BLACS_EXIT blacs_exit
#define BLACS_GRIDINFO blacs_gridinfo
#define BLACS_PNUM blacs_pnum
#define IGESD2D igesd2d
#define IGERV2D igerv2d
#define SGERV2D sgerv2d
#define SGESD2D sgesd2d
#define IGAMN2D igamn2d
#define SGSUM2D dgsum2d

#endif

/* Initialization functions */
#define my_Cblacs_pinfo Cblacs_pinfo
#define my_Cblacs_get Cblacs_get
#define my_Cblacs_gridinit Cblacs_gridinit
#define my_Cblacs_gridmap Cblacs_gridmap

/* Destruction functions */
#define my_Cblacs_gridexit Cblacs_gridexit
#define my_Cblacs_abort Cblacs_abort
#define my_Cblacs_exit Cblacs_exit

/* Information and miscellaneous */
#define my_Cblacs_gridinfo Cblacs_gridinfo
#define my_Cblacs_pnum Cblacs_pnum
#define my_Cblacs_barrier Cblacs_barrier

#define my_Cdgesd2d Cdgesd2d
#define my_Cdgerv2d Cdgerv2d

#ifdef FBLACS

extern "C"
{

    void my_Cblacs_pinfo(int* mypnum, int* nprocs);
    void my_Cblacs_get(int, int what, int* val);
    void my_Cblacs_barrier(int icontxt, char* scope);
    void my_Cblacs_gridinfo(
        int, int* nprow, int* npcol, int* myprow, int* mypcol);
    void my_Cblacs_gridinit(int* icontxt, char order[], int nprow, int npcol);
    void my_Cblacs_gridmap(
        int* icontxt, int* pmap, int ldpmap, int nprow, int npcol);
    void my_Cblacs_abort(int icontxt, int errornum);
    void my_Cblacs_gridexit(int icontxt);
    int my_Cblacs_pnum(int icontxt, int prow, int pcol);

    void my_Cdgesd2d(int, int, int, double*, int, int, int);
    void my_Cdgerv2d(int, int, int, double*, int, int, int);
    void Cigesd2d(int, int, int, int*, int, int, int);
}
#else

void Cigesd2d(int, int, int, int*, int, int, int);

extern "C"
{

    void Cblacs_pinfo(int* mypnum, int* nprocs);
    void Cblacs_get(int icontxt, int what, int* val);
    void Cblacs_barrier(int icontxt, char* scope);
    void Cblacs_gridinfo(
        int icontxt, int* nprow, int* npcol, int* myprow, int* mypcol);
    void Cblacs_gridinit(int* icontxt, char order[], int nprow, int npcol);
    void Cblacs_gridmap(
        int* icontxt, int* pmap, int ldpmap, int nprow, int npcol);
    void Cblacs_abort(int icontxt, int errornum);
    void Cblacs_gridexit(int icontxt);
    void Cblacs_exit(int idone);
    int Cblacs_pnum(int icontxt, int prow, int pcol);
    void Cdgesd2d(int, int, int, double*, int, int, int);
    void Cdgerv2d(int, int, int, double*, int, int, int);
    void Cdgsum2d(int, char*, char*, int, int, double*, int, int, int);
}

#endif

void my_Cdgsum2d(int, int, int, double*, int);

extern "C"
{

    void IGESD2D(int*, int*, int*, int*, int*, int*, int*);
    void SGESD2D(int*, int*, int*, double*, int*, int*, int*);
    void IGERV2D(int*, int*, int*, int*, int*, int*, int*);
    void SGERV2D(int*, int*, int*, double*, int*, int*, int*);
    void SGSUM2D(int*, CCHAR, CCHAR, int*, int*, double*, int*, int*, int*);

    void IGAMN2D(int*, CCHAR, CCHAR, int*, int*, int*, int*, int*, int*, int*,
        int*, int*);

    void BLACS_PINFO(int*, int*);
    void BLACS_GET(int*, int*, int*);
    void BLACS_BARRIER(int*, CCHAR);
    void BLACS_GRIDINFO(int*, int*, int*, int*, int*);
    void BLACS_GRIDINIT(int*, CCHAR, int*, int*);
    void BLACS_GRIDMAP(int*, int*, int*, int*, int*);
    void BLACS_ABORT(int*, int*);
    void BLACS_GRIDEXIT(int*);
    void BLACS_EXIT(int&);
    int BLACS_PNUM(int*, int*, int*);
}

#endif
