// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "mgmol_config.h"
#ifdef MGMOL_HAS_LIBROM

#include "LocGridOrbitals.h"
#include "MGmol.h"

#include "librom.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

template <class OrbitalsType>
const CAROM::Matrix* MGmol<OrbitalsType>::orbitals_to_carom_matrix(const OrbitalsType& orbitals)
{
    const int dim = orbitals.getLocNumpt();
    const int num_orbitals = orbitals.chromatic_number();

    CAROM::Options svd_options(dim, num_orbitals, 1);
    CAROM::BasisGenerator basis_generator(svd_options, false, "foo");

    for (int i = 0; i < num_orbitals; ++i)
        basis_generator.takeSample(orbitals.getPsi(i));

    return basis_generator.getSnapshotMatrix();
}

template <class OrbitalsType>
void MGmol<OrbitalsType>::carom_matrix_to_orbitals(const CAROM::Matrix* Psi, OrbitalsType& orbitals)
{
    Control& ct = *(Control::instance());
    Mesh* mesh = Mesh::instance();
    pb::GridFunc<ORBDTYPE> gf_psi(mesh->grid(), ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    CAROM::Vector psi;
    for (int i = 0; i < Psi->numColumns(); ++i)
    {
        Psi->getColumn(i, psi);
        gf_psi.assign(psi.getData());
        orbitals.setPsi(gf_psi, i);
    }
}

template <class OrbitalsType>
int MGmol<OrbitalsType>::save_orbital_snapshot(std::string file_path, OrbitalsType& orbitals)
{
    std::string snapshot_filename = file_path;
    struct stat s; 
    if (stat(file_path.c_str(), &s) == 0)
    {
        if (s.st_mode & S_IFDIR)
        {
            snapshot_filename = file_path + "/orbital";
        }
        else if (s.st_mode & S_IFREG)
        {
            snapshot_filename = file_path + "_orbital";
        }
        else
        {
            std::cout << file_path << " exists but is not a directory or a file." << std::endl;
            return 1;
        }
    }

    const int dim = orbitals.getLocNumpt();
    const int totalSamples = orbitals.chromatic_number();

    CAROM::Options svd_options(dim, totalSamples, 1);
    CAROM::BasisGenerator basis_generator(svd_options, false, snapshot_filename);

    for (int i = 0; i < totalSamples; ++i)
        basis_generator.takeSample(orbitals.getPsi(i));

    basis_generator.writeSnapshot();

    return 0;
}

template <class OrbitalsType>
void MGmol<OrbitalsType>::project_orbital(std::string file_path, int rdim, OrbitalsType& orbitals)
{
    const CAROM::Matrix* Psi = orbitals_to_carom_matrix(orbitals);

    CAROM::BasisReader reader(file_path);
    CAROM::Matrix* orbital_basis = reader.getSpatialBasis(rdim);

    CAROM::Matrix* Psi_reduced = orbital_basis->transposeMult(Psi);
    CAROM::Matrix* Psi_projected = orbital_basis->mult(Psi_reduced);
    
    carom_matrix_to_orbitals(Psi_projected, orbitals);
}

template class MGmol<LocGridOrbitals>;
template class MGmol<ExtendedGridOrbitals>;

#endif  // MGMOL_HAS_LIBROM
