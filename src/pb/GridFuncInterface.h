// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE


#ifndef GRIDFUNCINTERFACE_H_
#define GRIDFUNCINTERFACE_H_

#include "Timer.h"

namespace pb{

class GridFuncInterface
{
protected:
    static Timer   trade_bc_tm_;    
    static Timer   extend3D_tm_;    
    static Timer   restrict3D_tm_;    
    static Timer   prod_tm_;    
    static Timer   gather_tm_;    
    static Timer   scatter_tm_;    
    static Timer   all_gather_tm_;     
    static Timer  finishExchangeNorthSouth_tm_;
    static Timer  finishExchangeUpDown_tm_;
    static Timer  finishExchangeEastWest_tm_;

public:
    virtual ~GridFuncInterface(){}

    static void printTimers(std::ostream& os)
    {
       trade_bc_tm_.print(os);
       finishExchangeNorthSouth_tm_.print(os);
       finishExchangeUpDown_tm_.print(os);
       finishExchangeEastWest_tm_.print(os);
       extend3D_tm_.print(os);  
       restrict3D_tm_.print(os);
       prod_tm_.print(os);
       gather_tm_.print(os);
       scatter_tm_.print(os);
       all_gather_tm_.print(os); 
    } 

 
};

} // namespace pb
#endif
