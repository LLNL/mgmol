#ifndef FDkernels_H
#define FDkernels_H

#include "Grid.h"
#include "memory_space.h"

#include <iostream>

namespace pb
{

void printFDkernelTimers(std::ostream& os);

template <typename ScalarType>
void FDkernelDel2_2nd(const Grid& grid, ScalarType* v, ScalarType* b,
    const size_t nfunc, MemorySpace::Host);

template <typename ScalarType>
void FDkernelDel2_4th(const Grid& grid, ScalarType* v, ScalarType* b,
    const size_t nfunc, MemorySpace::Host);

template <typename ScalarType>
void FDkernelDel2_4th_Mehr(const Grid& grid, ScalarType* v, ScalarType* b,
    const size_t nfunc, MemorySpace::Host);

template <typename ScalarType>
void FDkernelDel2_6th(const Grid& grid, ScalarType* v, ScalarType* b,
    const size_t nfunc, MemorySpace::Host);

template <typename ScalarType>
void FDkernelDel2_8th(const Grid& grid, ScalarType* v, ScalarType* b,
    const size_t nfunc, MemorySpace::Host);

#ifdef HAVE_MAGMA
template <typename ScalarType>
void FDkernelDel2_2nd(const Grid& grid, ScalarType* v, ScalarType* b,
    const size_t nfunc, MemorySpace::Device)
{/*to be implemented*/ 
std::cerr<<"Function not implemented!!!"<<std::endl; abort();};

template <typename ScalarType>
void FDkernelDel2_4th(const Grid& grid, ScalarType* v, ScalarType* b,
    const size_t nfunc, MemorySpace::Device);

template <typename ScalarType>
void FDkernelDel2_4th_Mehr(const Grid& grid, ScalarType* v, ScalarType* b,
    const size_t nfunc, MemorySpace::Device)
{/*to be implemented*/
std::cerr<<"Function not implemented!!!"<<std::endl; abort();
};

template <typename ScalarType>
void FDkernelDel2_6th(const Grid& grid, ScalarType* v, ScalarType* b,
    const size_t nfunc, MemorySpace::Device)
{/*to be implemented*/
std::cerr<<"Function not implemented!!!"<<std::endl; abort();
};

template <typename ScalarType>
void FDkernelDel2_8th(const Grid& grid, ScalarType* v, ScalarType* b,
    const size_t nfunc, MemorySpace::Device)
{/*to be implemented*/
std::cerr<<"Function not implemented!!!"<<std::endl; abort();
};
#endif
}

#endif
