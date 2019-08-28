// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
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

using namespace std;

template <class T>
OrbitalsExtrapolation<T>::~OrbitalsExtrapolation()
{
    clearOldOrbitals();
}

template <class T>
void OrbitalsExtrapolation<T>::clearOldOrbitals()
{
    if (orbitals_minus1_ != nullptr)
    {
        delete orbitals_minus1_;
        orbitals_minus1_ = nullptr;
    }
#if EXTRAPOLATE_H
    delete hextrapol_;
#endif
}

template <class T>
bool OrbitalsExtrapolation<T>::getRestartData(T& orbitals)
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

template <class T>
void OrbitalsExtrapolation<T>::setupPreviousOrbitals(T** orbitals,
    ProjectedMatricesInterface* proj_matrices, LocalizationRegions* lrs,
    ClusterOrbitals* local_cluster, MasksSet* currentMasks, MasksSet* corrMasks,
    HDFrestart& h5f_file)
{
    if (onpe0)
        cout << "OrbitalsExtrapolation<T>::setupPreviousOrbitals()..." << endl;

    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    T* new_orbitals
        = new T("ForExtraploation", mygrid, mymesh->subdivx(), ct.numst, ct.bc,
            proj_matrices, lrs, currentMasks, corrMasks, local_cluster);

    new_orbitals->read_func_hdf5(h5f_file, "ExtrapolatedFunction");
    new_orbitals->incrementIterativeIndex();
    new_orbitals->incrementIterativeIndex();

    // swap pointers
    orbitals_minus1_ = (*orbitals);
    *orbitals        = new_orbitals;
}

template class OrbitalsExtrapolation<LocGridOrbitals>;
template class OrbitalsExtrapolation<ExtendedGridOrbitals>;
