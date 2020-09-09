// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "MVP_DMStrategy.h"
#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "MVPSolver.h"
#include "ProjectedMatricesInterface.h"

#include <vector>
using namespace std;

template <class OrbitalsType>
MVP_DMStrategy<OrbitalsType>::MVP_DMStrategy(MPI_Comm comm, ostream& os,
    Ions& ions, Rho<OrbitalsType>* rho, Energy<OrbitalsType>* energy,
    Electrostatic* electrostat, MGmol<OrbitalsType>* mgmol_strategy,
    OrbitalsType* orbitals, ProjectedMatricesInterface* proj_matrices,
    const bool use_old_dm)
    : orbitals_(orbitals),
      proj_matrices_(proj_matrices),
      comm_(comm),
      os_(os),
      ions_(ions),
      rho_(rho),
      energy_(energy),
      electrostat_(electrostat),
      global_indexes_(orbitals->getOverlappingGids()),
      mgmol_strategy_(mgmol_strategy),
      use_old_dm_(use_old_dm)
{
    assert(electrostat_ != nullptr);
    assert(energy_ != nullptr);
}

template <class OrbitalsType>
int MVP_DMStrategy<OrbitalsType>::update()
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout) << "MVP_DMStrategy<OrbitalsType>::update()..." << endl;
    }

    MVPSolver<OrbitalsType, dist_matrix::DistMatrix<DISTMATDTYPE>> solver(comm_,
        os_, ions_, rho_, energy_, electrostat_, mgmol_strategy_, ct.numst,
        ct.occ_width, global_indexes_, ct.dm_inner_steps, use_old_dm_);

    return solver.solve(*orbitals_);
}

template <class OrbitalsType>
void MVP_DMStrategy<OrbitalsType>::stripDM()
{
    if (use_old_dm_) proj_matrices_->stripDM();
}

template <class OrbitalsType>
void MVP_DMStrategy<OrbitalsType>::dressDM()
{
    if (use_old_dm_) proj_matrices_->dressupDM();
}

template class MVP_DMStrategy<LocGridOrbitals>;
template class MVP_DMStrategy<ExtendedGridOrbitals>;
