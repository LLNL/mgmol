// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "OrbitalsExtrapolation.h"
#include "ClusterOrbitals.h"
#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "MPIdata.h"
#include "MasksSet.h"
#include "Mesh.h"
#include "ProjectedMatricesInterface.h"

OrbitalsExtrapolation::~OrbitalsExtrapolation() { clearOldOrbitals(); }

void OrbitalsExtrapolation::clearOldOrbitals()
{
    if (orbitals_minus1_ != 0)
    {
        delete orbitals_minus1_;
        orbitals_minus1_ = 0;
    }
}

bool OrbitalsExtrapolation::getRestartData(LocGridOrbitals& orbitals)
{
    if (orbitals_minus1_ != 0)
    {
        orbitals.assign(*orbitals_minus1_);
        return true;
    }
    else
    {
        return false;
    }
}

void OrbitalsExtrapolation::setupPreviousOrbitals(LocGridOrbitals** orbitals,
    ProjectedMatricesInterface* proj_matrices, LocalizationRegions* lrs,
    ClusterOrbitals* local_cluster, MasksSet* currentMasks, MasksSet* corrMasks,
    HDFrestart& h5f_file)
{
    if (onpe0)
        cout << "OrbitalsExtrapolation::setupPreviousOrbitals()..." << endl;

    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    LocGridOrbitals* new_orbitals
        = new LocGridOrbitals(mygrid, mymesh->subdivx(), ct.numst, ct.bc,
            proj_matrices, lrs, currentMasks, corrMasks, local_cluster);

    new_orbitals->read_func_hdf5(h5f_file, "ExtrapolatedFunction");
    new_orbitals->incrementIterativeIndex();
    new_orbitals->incrementIterativeIndex();

    // swap pointers
    orbitals_minus1_ = (*orbitals);
    *orbitals        = new_orbitals;
}
