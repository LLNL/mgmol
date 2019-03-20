// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "NonOrthoDMStrategy.h"
#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "ProjectedMatricesInterface.h"

using namespace std;

template <class T>
NonOrthoDMStrategy<T>::NonOrthoDMStrategy(T* orbitals,
    ProjectedMatricesInterface* proj_matrices, const double mix)
    : orbitals_(orbitals), proj_matrices_(proj_matrices), mix_(mix)
{
}

template <class T>
void NonOrthoDMStrategy<T>::initialize()
{
    Control& ct = *(Control::instance());

    if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout) << "NonOrthoDMStrategy<T>::initialize()..." << endl;
    }
    proj_matrices_->updateDM(orbitals_->getIterativeIndex());
}

template <class T>
int NonOrthoDMStrategy<T>::update()
{
    assert(proj_matrices_ != 0);

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout) << "NonOrthoDMStrategy<T>::update() with mixing = "
                         << mix_ << endl;
    }

    // save old density matrix
    if (mix_ < 1.) proj_matrices_->saveDM();

    // compute new density matrix
    proj_matrices_->updateDM(orbitals_->getIterativeIndex());

    if (mix_ < 1.)
    {
        proj_matrices_->updateDMwithRelax(mix_, orbitals_->getIterativeIndex());
    }

    if (ct.verbose > 2)
    {
        double dd = proj_matrices_->getNel();
        if (onpe0)
            (*MPIdata::sout)
                << setprecision(8)
                << "test NonOrthoDMStrategy<T>::update(): Nel = " << dd << endl;
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

template class NonOrthoDMStrategy<LocGridOrbitals>;
template class NonOrthoDMStrategy<ExtendedGridOrbitals>;
