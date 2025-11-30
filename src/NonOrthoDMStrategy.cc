// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "NonOrthoDMStrategy.h"
#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "ProjectedMatricesInterface.h"

template <class OrbitalsType>
NonOrthoDMStrategy<OrbitalsType>::NonOrthoDMStrategy(
    ProjectedMatricesInterface* proj_matrices, const double mix)
    : proj_matrices_(proj_matrices), mix_(mix)
{
}

template <class OrbitalsType>
void NonOrthoDMStrategy<OrbitalsType>::initialize(OrbitalsType& /*orbitals*/)
{
    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    if (mmpi.PE0() && ct.verbose > 2)
    {
        (*MPIdata::sout) << "NonOrthoDMStrategy<T>::initialize()..."
                         << std::endl;
    }
    proj_matrices_->updateDM();
}

template <class OrbitalsType>
int NonOrthoDMStrategy<OrbitalsType>::update(OrbitalsType& /*orbitals*/)
{
    assert(proj_matrices_ != nullptr);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Control& ct     = *(Control::instance());

    if (mix_ <= 0.)
    {
        std::cerr << "NonOrthoDMStrategy, Invalid mixing value: " << mix_
                  << std::endl;
        MPI_Abort(mmpi.commSameSpin(), EXIT_FAILURE);
    }

    if (mmpi.PE0() && ct.verbose > 2)
    {
        (*MPIdata::sout) << "NonOrthoDMStrategy<T>::update() with mixing = "
                         << mix_ << std::endl;
    }

    // save old density matrix
    if (mix_ < 1.) proj_matrices_->saveDM();

    // compute new density matrix
    proj_matrices_->updateDM();

    if (mix_ < 1.)
    {
        proj_matrices_->updateDMwithRelax(mix_);
    }

    if (ct.verbose > 2)
    {
        double dd = proj_matrices_->getNel();
        if (mmpi.PE0())
            (*MPIdata::sout)
                << std::setprecision(8)
                << "test NonOrthoDMStrategy<T>::update(): Nel = " << dd
                << std::endl;
    }

    return 0; // success
}

template <class T>
void NonOrthoDMStrategy<T>::stripDM()
{
    if (mix_ < 1.) proj_matrices_->stripDM();
}

template <class T>
void NonOrthoDMStrategy<T>::dressDM()
{
    if (mix_ < 1.) proj_matrices_->dressupDM();
}

template class NonOrthoDMStrategy<LocGridOrbitals<ORBDTYPE>>;
template class NonOrthoDMStrategy<ExtendedGridOrbitals<ORBDTYPE>>;
