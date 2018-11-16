// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef RADIAL_FUNC_H
#define RADIAL_FUNC_H

namespace pb
{

double nu2(const double, const double, const double);
double nu4(double, double, double);
double nu6(double, double, double);
double epsilon6(const double, const double, const double);
double epsilon4(const double, const double);
double epsilon4(const double, const double, const double);
double epsilon2(const double, const double);
double epsilon2(const double, const double, const double);
double rho4(const double, const double);
double rho2(const double, const double);
double comp_potential(const double, const double);
double comp_charge(const double, const double);
double gaussian(const double r, const double rc);
double potential_of_gaussian(const double r, const double rc);
}

#endif
