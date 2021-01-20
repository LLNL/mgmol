#ifndef MGkernels_H
#define MGkernels_H

#include "Grid.h"

#include <iostream>

namespace pb
{

void printMGkernelTimers(std::ostream& os);

template <typename ScalarType>
void MGkernelExtend3D(ScalarType* coarse_data, const Grid& coarse_grid,
    ScalarType* fine_data, const Grid& fine_grid, const int nfunc);

template <typename ScalarType>
void MGkernelRestrict3D(ScalarType* fine_data, const Grid& fine_grid,
    ScalarType* coarse_data, const Grid& coarse_grid, const int nfunc);
}

#endif
