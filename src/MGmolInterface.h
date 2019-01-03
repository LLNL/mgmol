// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOLINTERFACE_H
#define MGMOLINTERFACE_H

class MGmolInterface
{
public:
    MGmolInterface(){}

    ~MGmolInterface(){}

    virtual int setupFromInput(const string input_file)=0;
    virtual int setupLRsFromInput(const string input_file)=0;
    virtual int setupConstraintsFromInput(const string input_file)=0;
    virtual void run()=0;
};

#endif

