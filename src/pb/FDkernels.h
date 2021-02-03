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
}

#endif
