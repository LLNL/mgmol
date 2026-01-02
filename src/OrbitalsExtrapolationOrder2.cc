// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "OrbitalsExtrapolationOrder2.h"
#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "ProjectedMatrices.h"

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrixTools.h"
#endif

template <class OrbitalsType>
void OrbitalsExtrapolationOrder2<OrbitalsType>::extrapolate_orbitals(
    OrbitalsType** orbitals, OrbitalsType* new_orbitals)
{
    Control& ct = *(Control::instance());

#ifdef MGMOL_USE_SCALAPACK
    bool use_dense_proj_mat = false;
    if (ct.OuterSolver() != OuterSolverType::ABPG
        && ct.OuterSolver() != OuterSolverType::NLCG)
    {
        ProjectedMatricesInterface* proj_matrices
            = (*orbitals)->getProjMatrices();
        if (dynamic_cast<
                ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>*>(
                proj_matrices))
            use_dense_proj_mat = true;
    }
#endif

    new_orbitals->assign(**orbitals);

    // do the extrapolation if previous orbitals exist (not at first step)
    if (OrbitalsExtrapolation<OrbitalsType>::orbitals_minus1_ != nullptr)
    {
        OrbitalsType* orbitals_minus1
            = OrbitalsExtrapolation<OrbitalsType>::orbitals_minus1_;
        if (ct.verbose > 1 && onpe0)
            (*MPIdata::sout) << "Extrapolate orbitals order 2..." << std::endl;

            // align orbitals_minus1_ with new_orbitals
#ifdef MGMOL_USE_SCALAPACK
        if (use_dense_proj_mat)
        {
            dist_matrix::DistMatrix<DISTMATDTYPE> matQ("Q", ct.numst, ct.numst);
            orbitals_minus1->computeGram(*new_orbitals, matQ);

            dist_matrix::DistMatrix<DISTMATDTYPE> yyt(
                "yyt", ct.numst, ct.numst);
            getProcrustesTransform(matQ, yyt);

            orbitals_minus1->multiply_by_matrix(matQ);
            orbitals_minus1->axpy((ORBDTYPE)-1., *new_orbitals);
            orbitals_minus1->multiply_by_matrix(yyt);
        }
        else // !use_dense_proj_mat
#endif
        {
            new_orbitals->scal(2.);
        }
        new_orbitals->axpy((ORBDTYPE)-1., *orbitals_minus1);

        delete orbitals_minus1;
    }

    // save data for next extrapolation
    if (ct.verbose > 1 && onpe0)
        (*MPIdata::sout) << "Set orbitals_minus1_ to values of orbitals"
                         << std::endl;
    OrbitalsExtrapolation<OrbitalsType>::orbitals_minus1_ = *orbitals;

    *orbitals = new_orbitals;

    (*orbitals)->incrementIterativeIndex();

    if (ct.isLocMode())
    {
        (*orbitals)->normalize();
        (*orbitals)->applyMask();
    }
    else
    {
        // DM (if not recomputed from scratch)
        // is consistant with orthonormal set of orbitals...
        if (ct.fullyOccupied()) (*orbitals)->orthonormalizeLoewdin();
    }
}

template class OrbitalsExtrapolationOrder2<LocGridOrbitals<ORBDTYPE>>;
template class OrbitalsExtrapolationOrder2<ExtendedGridOrbitals<ORBDTYPE>>;
