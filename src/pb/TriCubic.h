// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: TriCubic.h,v 1.2 2010/01/28 22:56:47 jeanluc Exp $
#ifndef PB_TRICUBIC_H
#define PB_TRICUBIC_H

#include "GridFunc.h"
#include <vector>

#ifdef HAVE_TRICUBIC

#include <mpi.h>

namespace pb
{

class Grid;
// class GridFunc;
template <class T>
class TriCubic
{
private:
    vector<vector<double>> spline_coeff_;

    const Grid& grid_;
    double h_[3];
    short bc_[3];

public:
    TriCubic(const Grid& grid, const short bc[3]);

    void computeSplineCoeffs(const T* const f);
    void computeSplineCoeffs(const GridFunc<T>& f, const GridFunc<T>& fx,
        const GridFunc<T>& fy, const GridFunc<T>& fz, const GridFunc<T>& fxy,
        const GridFunc<T>& fxz, const GridFunc<T>& fyz,
        const GridFunc<T>& fxyz);
    void getGradient(const double r[3], double dfdr[3], MPI_Comm);
    void getValues(const std::vector<double>&, std::vector<double>&, MPI_Comm);
};
}
#endif

#endif
