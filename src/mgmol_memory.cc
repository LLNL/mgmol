// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <stdio.h>
#include <stdlib.h>

#include "MPIdata.h"
#include "mgmol_memory.h"

#undef new

#define THRESHOLD4PRINT 10000
#define MAXNUMALLOCATIONS 100000
namespace MGmolMem
{
std::string Stamp::filename_ = "";
int Stamp::line_             = 0;
}

typedef struct
{
    long address;
    long size;
} Alloc_info;

Alloc_info* allocations[MAXNUMALLOCATIONS];
int nPos       = 0;
int total_size = 0;

void addTrack(long addr, long asize)
{
    static int high_water_mark = 0;

    Alloc_info* info  = (Alloc_info*)malloc(sizeof(Alloc_info));
    info->address     = addr;
    info->size        = asize;
    allocations[nPos] = info;
    nPos++;
    total_size += asize;
    if (total_size > high_water_mark)
    {
        high_water_mark = total_size;
        if (onpe0)
            printf(
                "MGmol allocation: new high_water_mark=%ld\n", high_water_mark);
    }

    if (asize > THRESHOLD4PRINT)
    {
        if (onpe0)
            printf("MGmol nPos %d: addr=%ld, size=%ld, filename=%s, line=%d\n",
                nPos, addr, asize, MGmolMem::Stamp::filename_.c_str(),
                MGmolMem::Stamp::line_);
    }

    if (nPos > MAXNUMALLOCATIONS)
    {
        printf("ERROR: Not enough memory slots!!!");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    // if(onpe0)printf("nPos=%d, addr=%ld, size=%ld\n",nPos,addr,asize);
}

bool removeTracking(long addr)
{
    bool bFound = false;
    long msize  = 0;

    for (int i = 0; i != nPos; i++)
    {
        if (allocations[i]->address == addr)
        {
            msize = allocations[i]->size;
            total_size -= msize;
            free(allocations[i]);
            bFound = true;
            for (int j = i; j < nPos - 1; ++j)
                allocations[j] = allocations[j + 1];
            nPos--;
            break;
        }
    }

    if (bFound && msize > THRESHOLD4PRINT)
        if (onpe0) printf("delete[%ld] of size %ld\n", addr, msize);

    if (!bFound)
        if (onpe0)
            printf(
                "WARNING: Unable to find allocation for delete[%ld]\n", addr);

    return bFound;
}

void* operator new(size_t sz) throw(std::bad_alloc)
{
    void* mem = malloc(sz);
    addTrack((long)mem, sz);
    if (mem)
    {
        return mem;
    }
    else
        throw std::bad_alloc();
}
void* operator new[](size_t sz) throw(std::bad_alloc)
{
    void* mem = malloc(sz);
    addTrack((long)mem, sz);
    if (mem == NULL) throw std::bad_alloc();
    return mem;
}

void operator delete(void* ptr) throw()
{
    removeTracking((long)ptr);
    free(ptr);
}
void operator delete[](void* ptr) throw()
{
    removeTracking((long)ptr);
    free(ptr);
}
