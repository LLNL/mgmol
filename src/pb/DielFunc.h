// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef DIEL_H
#define DIEL_H

// Dielfunc objects are essentially GridFunc<T> with periodic BC and members
// epsilon_max_, rho0_, drho0_
// They also have member functions to build the dielectric function
#include "GridFunc.h"

#include <cmath>

namespace pb
{
template <class T>
class DielFunc : public GridFunc<T>
{
private:
    static const T e0def_;
    static const T rho0def_;
    static const T drho0def_;

    T epsilon_max_;
    T rho0_;
    T drho0_;

public:
    // Constructor
    DielFunc(const Grid& my_grid) : GridFunc<T>(my_grid, 1, 1, 1)
    {
        epsilon_max_ = e0def_;
        rho0_        = rho0def_;
        drho0_       = drho0def_;
        // cout<<"Constructor for DielFunc at level "<<grid_.level()<<endl;
    };

    DielFunc(const Grid& my_grid, const T emax) : GridFunc<T>(my_grid, 1, 1, 1)
    {
        epsilon_max_ = emax;
        rho0_        = rho0def_;
        drho0_       = drho0def_;
        // cout<<"Constructor for DielFunc at level "<<grid_.level()<<endl;
    };

    DielFunc(const Grid& my_grid, const T emax, const T r0, const T dr0)
        : GridFunc<T>(my_grid, 1, 1, 1)
    {
        epsilon_max_ = emax;
        rho0_        = r0;
        drho0_       = dr0;
    };

    // copy constructor
    DielFunc(const DielFunc& A) : GridFunc<T>(A)
    {
        GridFunc<T>::bc_[0] = 1;
        GridFunc<T>::bc_[1] = 1;
        GridFunc<T>::bc_[2] = 1;
        epsilon_max_        = A.epsilon_max_;
        rho0_               = A.rho0_;
        drho0_              = A.drho0_;
    }

    // constructor
    DielFunc(const GridFunc<T>& A) : GridFunc<T>(A)
    {
        GridFunc<T>::bc_[0] = 1;
        GridFunc<T>::bc_[1] = 1;
        GridFunc<T>::bc_[2] = 1;
        epsilon_max_        = e0def_;
        rho0_               = rho0def_;
        drho0_              = drho0def_;
    }

    // copy constructor
    DielFunc(const DielFunc& A, Grid& my_grid) : GridFunc<T>(A, my_grid)
    {
        GridFunc<T>::bc_[0] = 1;
        GridFunc<T>::bc_[1] = 1;
        GridFunc<T>::bc_[2] = 1;
        epsilon_max_        = A.epsilon_max_;
        rho0_               = A.rho0_;
        drho0_              = A.drho0_;
    }

    ~DielFunc() {};

    DielFunc& operator=(const T);

    DielFunc& operator=(const DielFunc& A)
    {
        GridFunc<T>::bc_[0] = 1;
        GridFunc<T>::bc_[1] = 1;
        GridFunc<T>::bc_[2] = 1;
        epsilon_max_        = A.epsilon_max_;
        rho0_               = A.rho0_;
        drho0_              = A.drho0_;
        return *this;
    }

    double epsilon_rho(const T, const T, const T, const T);
    double epsilon_rho(const T rho);
    double depsilon_rho(const T, const T, const T, const T);
    double depsilon_rho(const T rho);

    void Gepsilon_rho(GridFunc<T>&, const T, const T);
    void Gepsilon_rho(GridFunc<T>&);
    void Gdepsilon_rho(GridFunc<T>&, GridFunc<T>&, const T, const T);
    void Gdepsilon_rho(GridFunc<T>&, GridFunc<T>&);

    double epsilon_max(void) const { return (double)epsilon_max_; }
};

template <class T>
const T DielFunc<T>::e0def_ = 78.36;
template <class T>
const T DielFunc<T>::rho0def_ = 0.0004;
template <class T>
const T DielFunc<T>::drho0def_ = 1.3;

} // namespace pb

#endif
