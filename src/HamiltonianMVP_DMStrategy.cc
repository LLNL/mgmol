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
#include "DistMatrix.h"
#include "HamiltonianMVPSolver.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesSparse.h"

template <class T1, class T2, class T3>
HamiltonianMVP_DMStrategy<T1, T2, T3>::HamiltonianMVP_DMStrategy(MPI_Comm comm,
    std::ostream& os, Ions& ions, Rho<T3>* rho, Energy<T3>* energy,
    Electrostatic* electrostat, MGmol<T3>* mgmol_strategy, T3* orbitals)
    : orbitals_(orbitals),
      comm_(comm),
      os_(os),
      ions_(ions),
      rho_(rho),
      energy_(energy),
      electrostat_(electrostat),
      global_indexes_(orbitals->getOverlappingGids()),
      mgmol_strategy_(mgmol_strategy)
{
    Control& ct = *(Control::instance());

    assert(electrostat_ != nullptr);
    assert(energy_ != nullptr);
    assert(ct.dm_inner_steps > 0);

    T2* projmatrices = dynamic_cast<T2*>(orbitals->getProjMatrices());

    solver_ = new HamiltonianMVPSolver<T1, T2, T3>(os_, ions_, rho_, energy_,
        electrostat_, mgmol_strategy_, ct.numst, ct.dm_inner_steps,
        projmatrices->getH(), true);
}

template <class T1, class T2, class T3>
HamiltonianMVP_DMStrategy<T1, T2, T3>::~HamiltonianMVP_DMStrategy()
{
    delete solver_;
}

template <class T1, class T2, class T3>
void HamiltonianMVP_DMStrategy<T1, T2, T3>::initialize()
{
}

template <class T1, class T2, class T3>
int HamiltonianMVP_DMStrategy<T1, T2, T3>::update()
{
    assert(solver_ != nullptr);

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout) << "HamiltonianMVP_DMStrategy::update()..."
                         << std::endl;
    }

    return solver_->solve(*orbitals_);
}

template <class T1, class T2, class T3>
void HamiltonianMVP_DMStrategy<T1, T2, T3>::stripDM()
{
}

template <class T1, class T2, class T3>
void HamiltonianMVP_DMStrategy<T1, T2, T3>::dressDM()
{
}

template <class T1, class T2, class T3>
void HamiltonianMVP_DMStrategy<T1, T2, T3>::reset()
{
    solver_->reset();
}

template class HamiltonianMVP_DMStrategy<dist_matrix::DistMatrix<DISTMATDTYPE>,
    ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>, LocGridOrbitals>;

template class HamiltonianMVP_DMStrategy<VariableSizeMatrix<sparserow>,
    ProjectedMatricesSparse, LocGridOrbitals>;

template class HamiltonianMVP_DMStrategy<dist_matrix::DistMatrix<DISTMATDTYPE>,
    ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>,
    ExtendedGridOrbitals>;
