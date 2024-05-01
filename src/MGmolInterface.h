// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOLINTERFACE_H
#define MGMOLINTERFACE_H

#include <cstring>

class MGmolInterface
{
public:
    MGmolInterface() {}

    virtual ~MGmolInterface() {}

    virtual int setupFromInput(const std::string input_file)            = 0;
    virtual int setupLRs(const std::string input_file)                  = 0;
    virtual int setupConstraintsFromInput(const std::string input_file) = 0;
    virtual void setup()                                                = 0;
    virtual void run()                                                  = 0;
    
};

#endif
