// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef LAPFACTORY_H
#define LAPFACTORY_H

#include "Laph2.h"
#include "Laph4.h"
#include "Laph6.h"
#include "Laph8.h"
#include "Laph4M.h"
#include "Laph4MP.h"
#include "MPIdata.h"

template <class T>
class LapFactory
{
public:
    static pb::Lap<T>* createLap(const pb::Grid& grid, const int type)
    {
        pb::Lap<T>* lap;
        switch( type ){
            case 0:
                lap = new pb::Laph4M<T>(grid);
                break;
            case 1:
                lap = new pb::Laph2<T>(grid);
                break;
            case 2:
                lap = new pb::Laph4<T>(grid);
                break;
            case 3:
                lap = new pb::Laph6<T>(grid);
                break;
            case 4:
                lap = new pb::Laph8<T>(grid);
                break;
            case 10:
                lap = new pb::Laph4MP<T>(grid);
                break;
            default:
                (*MPIdata::serr)<<"LapFactory::createLap() --- option invalid:"
                    <<type<<std::endl;
                exit(2);
        }
        return lap;
    }
};

#endif
