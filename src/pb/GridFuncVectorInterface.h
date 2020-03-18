// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef GRIDFUNCVECTORINTERFACE_H_
#define GRIDFUNCVECTORINTERFACE_H_

#include "Timer.h"

namespace pb
{

class GridFuncVectorInterface
{
protected:
    static Timer trade_bc_tm_;
    static Timer trade_bc_colors_tm_;
    static Timer prod_tm_;
    static Timer finishExchangeNorthSouth_tm_;
    static Timer finishExchangeUpDown_tm_;
    static Timer finishExchangeEastWest_tm_;
    static Timer wait_north_south_tm_;
    static Timer wait_up_down_tm_;
    static Timer wait_east_west_tm_;

public:
    virtual ~GridFuncVectorInterface() {}

    static void printTimers(std::ostream& os)
    {
        trade_bc_tm_.print(os);
        trade_bc_colors_tm_.print(os);
        prod_tm_.print(os);
        wait_north_south_tm_.print(os);
        wait_up_down_tm_.print(os);
        wait_east_west_tm_.print(os);
        finishExchangeNorthSouth_tm_.print(os);
        finishExchangeUpDown_tm_.print(os);
        finishExchangeEastWest_tm_.print(os);
    }
};
}
#endif
