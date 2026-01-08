// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "EigenDMStrategy.h"
#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "ProjectedMatrices.h"

template <class OrbitalsType>
EigenDMStrategy<OrbitalsType>::EigenDMStrategy(
    ProjectedMatricesInterface* proj_matrices)
    : proj_matrices_(proj_matrices)
{
}

template <class OrbitalsType>
void EigenDMStrategy<OrbitalsType>::initialize(OrbitalsType& orbitals)
{
    update(orbitals);
}

template <class OrbitalsType>
int EigenDMStrategy<OrbitalsType>::update(OrbitalsType& orbitals)
{
#ifdef MGMOL_USE_SCALAPACK
    Control& ct = *(Control::instance());

    dist_matrix::DistMatrix<DISTMATDTYPE> zz("Z", ct.numst, ct.numst);

    ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>* pmat
        = dynamic_cast<
            ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>*>(
            proj_matrices_);
    pmat->updateDMwithEigenstatesAndRotate(zz);

    // if( onpe0 && ct.verbose>2 )
    //    (*MPIdata::sout)<<"get_dm_diag: rotate orbitals "<<endl;
    orbitals.multiply_by_matrix(zz);
    orbitals.setDataWithGhosts();
    orbitals.trade_boundaries();
#else
    (void)orbitals;
    std::cerr << "EigenDMStrategy<OrbitalsType>::update(OrbitalsType& "
                 "orbitals) not implemented"
              << std::endl;
    abort();
#endif

    return 0;
}

template class EigenDMStrategy<LocGridOrbitals<ORBDTYPE>>;
template class EigenDMStrategy<ExtendedGridOrbitals<ORBDTYPE>>;
