// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef MatricesBlacsContext_H
#define MatricesBlacsContext_H

#include "BlacsContext.h"
#include <cassert>
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif


class MatricesBlacsContext
{
    MatricesBlacsContext()
       :blactxt_(0),
        comm_(-1),
        size_(-1),
        n_(-1)
    {}
    
    ~MatricesBlacsContext()
    {
        if( blactxt_!=0 )delete blactxt_;
    };
    
    dist_matrix::BlacsContext *blactxt_;
    MPI_Comm comm_;
    int size_;
    int n_;

public:
    static MatricesBlacsContext& instance();
    
    void setup(const MPI_Comm comm, const int size);
    
    dist_matrix::BlacsContext* bcxt()const
    {
        assert( blactxt_!=0 );
        return blactxt_;
    }
    
    void print(std::ostream& os)const
    {
        os<<"MatricesBlacsContext: "<<n_<<" x "<<n_<<std::endl;
    }
    
    void clear()
    {
        if( blactxt_!=0 )
        {
            delete blactxt_;
            blactxt_=0;
        }
        comm_=-1;
        size_=-1;
        n_=-1;
    }
};

#endif
