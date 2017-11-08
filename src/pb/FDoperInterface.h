// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
 * FDoperInterface class for timing variables
*/
#ifndef _FDOPERINTERFACE_H_
#define _FDOPERINTERFACE_H_

#include "Timer.h"
#include <iostream>

namespace pb{

class FDoperInterface{
protected:
    static Timer   del2_4th_Mehr_tm_;    
    static Timer   rhs_4th_Mehr1_tm_;    
    static Timer   rhs_tm_;    
    static Timer   del2_2nd_tm_;    
    static Timer   del2_4th_tm_;  
    static Timer   del2_4th_wpot_tm_;  

public:
    virtual ~FDoperInterface(){}
    
    static void printTimers(std::ostream& os)
    {
       del2_4th_Mehr_tm_.print(os);    
       rhs_tm_.print(os);    
       rhs_4th_Mehr1_tm_.print(os);    
       del2_2nd_tm_.print(os);    
       del2_4th_tm_.print(os); 
       del2_4th_wpot_tm_.print(os); 
    }
};

} // namespace
#endif  