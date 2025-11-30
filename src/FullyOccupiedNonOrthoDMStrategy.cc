// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "FullyOccupiedNonOrthoDMStrategy.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "ProjectedMatricesInterface.h"

template <class OrbitalsType>
FullyOccupiedNonOrthoDMStrategy<OrbitalsType>::FullyOccupiedNonOrthoDMStrategy(
    ProjectedMatricesInterface* proj_matrices)
    : proj_matrices_(proj_matrices)
{
}

template <class OrbitalsType>
void FullyOccupiedNonOrthoDMStrategy<OrbitalsType>::initialize(
    OrbitalsType& orbitals)
{
    update(orbitals);
}

template <class OrbitalsType>
int FullyOccupiedNonOrthoDMStrategy<OrbitalsType>::update(
    OrbitalsType& orbitals)
{
    assert(proj_matrices_ != nullptr);
    (void)orbitals;

    proj_matrices_->setDMto2InvS();

#if 0 /* TEST */
    int dd=proj_matrices_->getNel();
    if( onpe0 && ct.verbose>2 ){
        (*MPIdata::sout)<<"test get_dm_diag:  Nel = "<<dd<<endl;
    }
    assert( (dd-nel)==0 );
#endif

    return 0; // success
}

template class FullyOccupiedNonOrthoDMStrategy<LocGridOrbitals<ORBDTYPE>>;
template class FullyOccupiedNonOrthoDMStrategy<ExtendedGridOrbitals<ORBDTYPE>>;
