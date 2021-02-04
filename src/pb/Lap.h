// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

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
public:
    using memory_space_type = typename FDoperInterface::memory_space_type;

protected:
    std::string name_;

public:
    Lap(const Grid& mygrid) : FDoper<T>(mygrid) {}

    Lap& operator=(const Lap& v) = default;

    // A->B
    void apply(GridFunc<T>& A, GridFunc<T>& B) override = 0;
    virtual void applyWithPot(GridFunc<T>&, const double* const, T*)
    {
        std::cerr << "ERROR: Lap::applyWithPot() not implemented" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    virtual void apply(GridFuncVector<T, memory_space_type>& A,
        GridFuncVector<T, memory_space_type>& B)
        = 0;

    std::string name() const { return name_; }

    virtual void jacobi(GridFunc<T>&, const GridFunc<T>&, GridFunc<T>&) = 0;
    virtual void jacobi(GridFuncVector<T, memory_space_type>&,
        const GridFuncVector<T, memory_space_type>&, GridFunc<T>&)
        = 0;
    virtual void jacobi(GridFuncVector<T, memory_space_type>&,
        const GridFuncVector<T, memory_space_type>&,
        GridFuncVector<T, memory_space_type>&)
        = 0;

    void jacobi(GridFunc<T>&, const GridFunc<T>&, GridFunc<T>&, const double);
    void jacobi(GridFuncVector<T, memory_space_type>&,
        const GridFuncVector<T, memory_space_type>&, GridFunc<T>&,
        const double);
    void jacobi(GridFuncVector<T, memory_space_type>&,
        const GridFuncVector<T, memory_space_type>&,
        GridFuncVector<T, memory_space_type>&, const double);

    double energyES(GridFunc<T>&, GridFunc<T>&);
    virtual double diagEl(void) const    = 0;
    virtual double invDiagEl(void) const = 0;
    virtual void setLowerOrderGrid(void) = 0;

    ~Lap() override {}
};

} // namespace pb

#endif
