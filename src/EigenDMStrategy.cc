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

template <class T>
EigenDMStrategy<T>::EigenDMStrategy(
    T* current_orbitals, ProjectedMatricesInterface* proj_matrices)
    : current_orbitals_(current_orbitals), proj_matrices_(proj_matrices)
{
}

template <class T>
void EigenDMStrategy<T>::initialize()
{
    update();
}

template <class T>
int EigenDMStrategy<T>::update()
{
    Control& ct = *(Control::instance());

    dist_matrix::DistMatrix<DISTMATDTYPE> zz("Z", ct.numst, ct.numst);

    ProjectedMatrices* pmat = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    pmat->updateDMwithEigenstatesAndRotate(
        current_orbitals_->getIterativeIndex(), zz);

    // if( onpe0 && ct.verbose>2 )
    //    (*MPIdata::sout)<<"get_dm_diag: rotate orbitals "<<endl;
    current_orbitals_->multiply_by_matrix(zz);
    current_orbitals_->setDataWithGhosts();
    current_orbitals_->trade_boundaries();

    return 0;
}

template class EigenDMStrategy<LocGridOrbitals>;
template class EigenDMStrategy<ExtendedGridOrbitals>;
