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
#include "ReplicatedMatrix.h"

#include <vector>
using namespace std;

template <class OrbitalsType, class MatrixType>
MVP_DMStrategy<OrbitalsType, MatrixType>::MVP_DMStrategy(MPI_Comm comm,
    ostream& os, Ions& ions, Rho<OrbitalsType>* rho,
    Energy<OrbitalsType>* energy, Electrostatic* electrostat,
    Hamiltonian<OrbitalsType>* hamiltonian, MGmol<OrbitalsType>* mgmol_strategy,
    const std::vector<std::vector<int>>& overlappingGids,
    ProjectedMatricesInterface* proj_matrices, const bool use_old_dm)
    : proj_matrices_(proj_matrices),
      comm_(comm),
      os_(os),
      ions_(ions),
      rho_(rho),
      energy_(energy),
      electrostat_(electrostat),
      hamiltonian_(hamiltonian),
      global_indexes_(overlappingGids),
      mgmol_strategy_(mgmol_strategy),
      use_old_dm_(use_old_dm)
{
    assert(electrostat_ != nullptr);
    assert(energy_ != nullptr);
}

template <class OrbitalsType, class MatrixType>
int MVP_DMStrategy<OrbitalsType, MatrixType>::update(OrbitalsType& orbitals)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout)
            << "MVP_DMStrategy<OrbitalsType,MatrixType>::update()..." << endl;
    }

    MVPSolver<OrbitalsType, MatrixType> solver(comm_, os_, ions_, rho_, energy_,
        electrostat_, hamiltonian_, mgmol_strategy_, ct.numst, ct.occ_width,
        global_indexes_, ct.dm_inner_steps, ct.dm_mix, ct.dm_tol, use_old_dm_);

    return solver.solve(orbitals);
}

template <class OrbitalsType, class MatrixType>
void MVP_DMStrategy<OrbitalsType, MatrixType>::stripDM()
{
    if (use_old_dm_) proj_matrices_->stripDM();
}

template <class OrbitalsType, class MatrixType>
void MVP_DMStrategy<OrbitalsType, MatrixType>::dressDM()
{
    if (use_old_dm_) proj_matrices_->dressupDM();
}

#ifdef MGMOL_USE_SCALAPACK
template class MVP_DMStrategy<LocGridOrbitals<ORBDTYPE>,
    dist_matrix::DistMatrix<double>>;
template class MVP_DMStrategy<ExtendedGridOrbitals<ORBDTYPE>,
    dist_matrix::DistMatrix<double>>;
#endif
template class MVP_DMStrategy<LocGridOrbitals<ORBDTYPE>, ReplicatedMatrix>;
template class MVP_DMStrategy<ExtendedGridOrbitals<ORBDTYPE>, ReplicatedMatrix>;
