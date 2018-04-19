// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "MatricesBlacsContext.h"
#include "DistMatrix.h"

using namespace std;

// allows zero or one object
MatricesBlacsContext& MatricesBlacsContext::instance()
{
    static MatricesBlacsContext instance;
    return instance;
}

void MatricesBlacsContext::setup(const MPI_Comm comm, const int size)
{
    assert( size>=0 );
    
    comm_=comm;
    size_=size;
    
    if( size_>=0 )
    {
        const int bsize=dist_matrix::DistMatrix<DISTMATDTYPE>::getBlockSize();
        n_=(int)ceilf((float)size_
                /(float)(3*bsize));
        if( size_>2*bsize )n_=max(n_,2);
        blactxt_= new dist_matrix::BlacsContext(comm_,'s',n_*n_);
    }
    else
    {
        blactxt_=0;
    }
};
