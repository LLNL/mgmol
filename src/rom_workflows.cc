// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "rom_workflows.h"

template <class OrbitalsType>
void readRestartFiles(MGmolInterface *mgmol_)
{
    Control& ct              = *(Control::instance());
    ROMPrivateOptions rom_options = ct.getROMOptions();

    MGmol<OrbitalsType> *mgmol = static_cast<MGmol<OrbitalsType> *>(mgmol_);

    OrbitalsType *orbitals = mgmol->loadOrbitalFromRestartFile(rom_options.restart_filename);

    mgmol->save_orbital_snapshot("test", *orbitals);

    delete orbitals;
}

template void readRestartFiles<LocGridOrbitals>(MGmolInterface *mgmol_);
template void readRestartFiles<ExtendedGridOrbitals>(MGmolInterface *mgmol_);