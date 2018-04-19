// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef RADIALINTER_H
#define RADIALINTER_H

#include "RadialFunction.h"

class RadialInter:public RadialFunction
{
private:
    

public:

    RadialInter(const std::vector<double>& x):RadialFunction(x)
    {
    }

    RadialInter():RadialFunction()
    {
    }

    RadialInter(std::vector<double>& x, std::vector<std::vector<double> >& y):
        RadialFunction(x,y)
    {
    }

    double linint(const double x, const int j=0)const;
    double cubint(const double x, const int j=0)const;
};

#endif
