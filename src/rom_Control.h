// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef ROM_CONTROL_H
#define ROM_CONTROL_H

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

enum class ROMStage
{
    OFFLINE,
    ONLINE,
    RESTORE,
    UNSUPPORTED
};

/* Stored as a private member variable of Control class */
struct ROMPrivateOptions
{
    std::string restart_filename = "";
    ROMStage rom_stage = ROMStage::UNSUPPORTED;
};

#endif  // ROM_CONTROL_H
