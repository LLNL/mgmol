// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "romOffline.h"
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"

#include <string>
#include <iostream>
#include <fstream>

// Save the wavefunction snapshots
template <class OrbitalsType>
int save_orbital_snapshot(std::string snapshot_dir, OrbitalsType& orbitals)
{
    std::string snapshot_filename = snapshot_dir + "/orbital_snapshot";
    std::cout << "Creating snapshot file " << snapshot_filename << std::endl;

    const int dim = orbitals.numst();
    const int totalSamples = orbitals.chromatic_number();
    std::cout << "Snapshot size: " << dim << "*" << totalSamples << std::endl;

    CAROM::Options svd_options(dim, totalSamples, 1);
    CAROM::BasisGenerator basis_generator(svd_options, false, snapshot_filename);

    for (int i = 0; i < totalSamples; ++i)
        basis_generator.takeSample(orbitals.psi(i));

    basis_generator.writeSnapshot();

    return 0;
}
