// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <iostream>
#include "MPIdata.h"

class Index
{
    int i_;
    int j_;
public:
    Index(const int i, const int j)
    {
        i_=i;
        j_=j;
    }
    
    int i()const{ return i_;}
    int j()const{ return j_;}
    
    void print()const
    {
        (*MPIdata::sout)<<"("<<i_<<","<<j_<<")"<<std::endl;
    }

};

bool operator<(const Index i, const Index j)
{
    bool alpha=(i.i()==j.i());
    if( alpha )return i.j()<j.j();
    else       return i.i()<j.i();
}

