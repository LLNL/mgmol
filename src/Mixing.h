// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$

#ifndef MIXING
#define MIXING


template <class T>
class Mixing
{
public:
    Mixing(){};
    
    virtual ~Mixing(){};
    
    virtual void update(T& res,T& work)=0;
    virtual void restart(void)=0;
};

#endif
