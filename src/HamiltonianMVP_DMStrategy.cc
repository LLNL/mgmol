// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "HamiltonianMVP_DMStrategy.h"
#include "Control.h"
#include "HamiltonianMVPSolver.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesSparse.h"
#include "ReplicatedMatrix.h"

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrix.h"
#endif

template <class MatrixType, class ProjMatrixType, class OrbitalsType>
HamiltonianMVP_DMStrategy<MatrixType, ProjMatrixType,
    OrbitalsType>::HamiltonianMVP_DMStrategy(MPI_Comm comm, std::ostream& os,
    Ions& ions, Rho<OrbitalsType>* rho, Energy<OrbitalsType>* energy,
    Electrostatic* electrostat, Hamiltonian<OrbitalsType>* hamiltonian,
    MGmol<OrbitalsType>* mgmol_strategy, OrbitalsType* orbitals)
    : comm_(comm),
      os_(os),
      ions_(ions),
      rho_(rho),
      energy_(energy),
      electrostat_(electrostat),
      hamiltonian_(hamiltonian),
      global_indexes_(orbitals->getOverlappingGids()),
      mgmol_strategy_(mgmol_strategy)
{
    Control& ct = *(Control::instance());

    assert(electrostat_ != nullptr);
    assert(energy_ != nullptr);
    assert(ct.dm_inner_steps > 0);

    ProjMatrixType* projmatrices
        = dynamic_cast<ProjMatrixType*>(orbitals->getProjMatrices());

    solver_
        = new HamiltonianMVPSolver<MatrixType, ProjMatrixType, OrbitalsType>(
            os_, ions_, rho_, energy_, electrostat_, hamiltonian_,
            mgmol_strategy_, ct.numst, ct.dm_inner_steps, projmatrices->getH(),
            true);
}

template <class MatrixType, class ProjMatrixType, class OrbitalsType>
HamiltonianMVP_DMStrategy<MatrixType, ProjMatrixType,
    OrbitalsType>::~HamiltonianMVP_DMStrategy()
{
    delete solver_;
}

template <class MatrixType, class ProjMatrixType, class OrbitalsType>
void HamiltonianMVP_DMStrategy<MatrixType, ProjMatrixType,
    OrbitalsType>::initialize(OrbitalsType&)
{
}

template <class MatrixType, class ProjMatrixType, class OrbitalsType>
int HamiltonianMVP_DMStrategy<MatrixType, ProjMatrixType, OrbitalsType>::update(
    OrbitalsType& orbitals)
{
    assert(solver_ != nullptr);

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout) << "HamiltonianMVP_DMStrategy::update()..."
                         << std::endl;
    }

    return solver_->solve(orbitals);
}

template <class MatrixType, class ProjMatrixType, class OrbitalsType>
void HamiltonianMVP_DMStrategy<MatrixType, ProjMatrixType,
    OrbitalsType>::stripDM()
{
}

template <class MatrixType, class ProjMatrixType, class OrbitalsType>
void HamiltonianMVP_DMStrategy<MatrixType, ProjMatrixType,
    OrbitalsType>::dressDM()
{
}

template <class MatrixType, class ProjMatrixType, class OrbitalsType>
void HamiltonianMVP_DMStrategy<MatrixType, ProjMatrixType,
    OrbitalsType>::reset()
{
    solver_->reset();
}

#ifdef MGMOL_USE_SCALAPACK
template class HamiltonianMVP_DMStrategy<dist_matrix::DistMatrix<DISTMATDTYPE>,
    ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>,
    LocGridOrbitals<ORBDTYPE>>;
template class HamiltonianMVP_DMStrategy<dist_matrix::DistMatrix<DISTMATDTYPE>,
    ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>,
    ExtendedGridOrbitals<ORBDTYPE>>;
#endif
template class HamiltonianMVP_DMStrategy<VariableSizeMatrix<sparserow>,
    ProjectedMatricesSparse, LocGridOrbitals<ORBDTYPE>>;
template class HamiltonianMVP_DMStrategy<ReplicatedMatrix,
    ProjectedMatrices<ReplicatedMatrix>, ExtendedGridOrbitals<ORBDTYPE>>;
