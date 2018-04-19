// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef DMSTRATEGY_H
#define DMSTRATEGY_H

class DMStrategy
{
public:
    virtual void initialize()=0;
    virtual int update()=0;
    
    virtual ~DMStrategy(){};
    
    // tells if strategy needs an updated H matrix
    // to update DM
    virtual bool needH()const=0;
    
    virtual void stripDM()=0;
    virtual void dressDM()=0;
    
    /*!
     * reset initial DM to the one used in previous update() call
     */
    virtual void reset()=0;
};

#endif
