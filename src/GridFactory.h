// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef GRIDFACTORY_H
#define GRIDFACTORY_H

#include "Grid.h"
#include "MPIdata.h"

class GridFactory
{
public:
    static pb::Grid* createGrid(const unsigned ngpts[3], const double origin[3],
        const double lattice[3], const int lap_type, const bool diel_flag,
        const pb::PEenv& myPEenv)
    {
        int nghosts = 0;
        switch (lap_type)
        {
            case 0:
                nghosts = diel_flag ? 2 : 1;
                break;
            case 1:
                nghosts = 1;
                break;
            case 2:
                nghosts = 2;
                break;
            case 3:
                nghosts = 3;
                break;
            case 4:
                nghosts = 4;
                break;
            case 10:
                nghosts = diel_flag ? 2 : 1;
                break;
            default:
                (*MPIdata::serr) << "lap_type = " << lap_type << std::endl;
                (*MPIdata::serr)
                    << "GridFactory::createGrid() --- option invalid."
                    << std::endl;
                exit(2);
        }
        return (new pb::Grid(origin, lattice, ngpts, myPEenv, nghosts));
    }
};

#endif
