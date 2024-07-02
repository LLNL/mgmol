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
    RESTORE,    // TODO(kevin): what stage is this?
    BUILD,
    UNSUPPORTED
};

enum class ROMVariable
{
    ORBITALS,
    POTENTIAL,
    NONE
};

/* Stored as a private member variable of Control class */
struct ROMPrivateOptions
{
    ROMStage rom_stage = ROMStage::UNSUPPORTED;

    std::string restart_file_fmt = "";
    int restart_file_minidx = -1;
    int restart_file_maxidx = -1;
    std::string basis_file = "";
    ROMVariable variable=ROMVariable::NONE;

    /* save librom snapshot matrix at FOM simulation. */
    bool save_librom_snapshot = false;

    /* options for ROM building */
    int num_potbasis = -1;
};

#endif  // ROM_CONTROL_H
