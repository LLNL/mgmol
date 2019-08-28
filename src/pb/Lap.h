// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Lap.h,v 1.11 2009/07/16 23:36:28 jeanluc Exp $
#ifndef PB_LAP_H
#define PB_LAP_H

#include "FDoper.h"
#include "GridFuncVector.h"

#include <string>

namespace pb
{
template <class T>
class Lap : public FDoper<T>
{

protected:
    std::string name_;

public:
    Lap(const Grid& mygrid) : FDoper<T>(mygrid) {}

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B) override = 0;
    virtual void applyWithPot(GridFunc<T>& A, const double* const, T*)
    {
        std::cerr << "ERROR: Lap::applyWithPot() not implemented" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    virtual void apply(GridFuncVector<T>& A, GridFuncVector<T>& B) = 0;

    std::string name() const { return name_; }

    virtual void jacobi(GridFunc<T>&, const GridFunc<T>&, GridFunc<T>&) = 0;
    virtual void jacobi(
        GridFuncVector<T>&, const GridFuncVector<T>&, GridFunc<T>&)
        = 0;
    virtual void jacobi(
        GridFuncVector<T>&, const GridFuncVector<T>&, GridFuncVector<T>&)
        = 0;

    void jacobi(GridFunc<T>&, const GridFunc<T>&, GridFunc<T>&, const double);
    void jacobi(GridFuncVector<T>&, const GridFuncVector<T>&, GridFunc<T>&,
        const double);
    void jacobi(GridFuncVector<T>&, const GridFuncVector<T>&,
        GridFuncVector<T>&, const double);

    double energyES(GridFunc<T>&, GridFunc<T>&);
    virtual double diagEl(void) const    = 0;
    virtual double invDiagEl(void) const = 0;
    virtual void setLowerOrderGrid(void) = 0;

    ~Lap() override {}
};

} // namespace pb

#endif
