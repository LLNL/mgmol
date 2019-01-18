// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "MVP_DMStrategy.h"
#include "Control.h"
#include "Ions.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "MVPSolver.h"
#include "ProjectedMatricesInterface.h"

#include <vector>
using namespace std;

template <class T>
MVP_DMStrategy<T>::MVP_DMStrategy(MPI_Comm comm, ostream& os, Ions& ions,
    Rho<T>* rho,
    Energy<T>* energy, Electrostatic* electrostat, MGmol<T>* mgmol_strategy,
    T* orbitals, ProjectedMatricesInterface* proj_matrices,
    const bool use_old_dm)
    : comm_(comm),
      os_(os),
      ions_(ions),
      rho_(rho),
      energy_(energy),
      electrostat_(electrostat),
      global_indexes_(orbitals->getOverlappingGids()),
      orbitals_(orbitals),
      proj_matrices_(proj_matrices),
      mgmol_strategy_(mgmol_strategy),
      use_old_dm_(use_old_dm)
{
    assert(electrostat_ != 0);
    assert(energy_ != 0);
}

template <class T>
int MVP_DMStrategy<T>::update()
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout) << "MVP_DMStrategy<T>::update()..." << endl;
    }

    MVPSolver<T> solver(comm_, os_, ions_, rho_, energy_,
        electrostat_,
        mgmol_strategy_, ct.numst, ct.occ_width, ct.getNel(), global_indexes_,
        ct.dm_inner_steps, use_old_dm_);

    return solver.solve(*orbitals_);
}

template <class T>
void MVP_DMStrategy<T>::stripDM()
{
    if (use_old_dm_) proj_matrices_->stripDM();
}

template <class T>
void MVP_DMStrategy<T>::dressDM()
{
    if (use_old_dm_) proj_matrices_->dressupDM();
}

template class MVP_DMStrategy<LocGridOrbitals>;
template class MVP_DMStrategy<ExtendedGridOrbitals>;
