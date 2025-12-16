// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "OrbitalsExtrapolation.h"
#include "ClusterOrbitals.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "MPIdata.h"
#include "MasksSet.h"
#include "Mesh.h"
#include "ProjectedMatricesInterface.h"

template <class OrbitalsType>
OrbitalsExtrapolation<OrbitalsType>::~OrbitalsExtrapolation()
{
    if (orbitals_minus1_ != nullptr)
    {
        delete orbitals_minus1_;
        orbitals_minus1_ = nullptr;
    }
}

template <class OrbitalsType>
void OrbitalsExtrapolation<OrbitalsType>::clearOldOrbitals()
{
    if (orbitals_minus1_ != nullptr)
    {
        delete orbitals_minus1_;
        orbitals_minus1_ = nullptr;
    }
}

template <class OrbitalsType>
bool OrbitalsExtrapolation<OrbitalsType>::getRestartData(OrbitalsType& orbitals)
{
    if (orbitals_minus1_ != nullptr)
    {
        orbitals.assign(*orbitals_minus1_);
        return true;
    }
    else
    {
        return false;
    }
}

template <class OrbitalsType>
void OrbitalsExtrapolation<OrbitalsType>::setupPreviousOrbitals(
    OrbitalsType** orbitals, ProjectedMatricesInterface* proj_matrices,
    std::shared_ptr<LocalizationRegions> lrs, ClusterOrbitals* local_cluster,
    MasksSet* currentMasks, MasksSet* corrMasks, HDFrestart& h5f_file)
{
    if (onpe0)
        std::cout
            << "OrbitalsExtrapolation<OrbitalsType>::setupPreviousOrbitals()..."
            << std::endl;

    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    OrbitalsType* new_orbitals = new OrbitalsType("ForExtraploation", mygrid,
        mymesh->subdivx(), ct.numst, ct.bcWF, proj_matrices, lrs, currentMasks,
        corrMasks, local_cluster);

    new_orbitals->read_func_hdf5(h5f_file, "ExtrapolatedFunction");
    new_orbitals->incrementIterativeIndex();
    new_orbitals->incrementIterativeIndex();

    // swap pointers
    orbitals_minus1_ = (*orbitals);
    *orbitals        = new_orbitals;
}

template class OrbitalsExtrapolation<LocGridOrbitals<ORBDTYPE>>;
template class OrbitalsExtrapolation<ExtendedGridOrbitals<ORBDTYPE>>;
