// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "OrbitalsExtrapolationOrder3.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "ProjectedMatrices.h"

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrixTools.h"
#endif

template <class OrbitalsType>
void OrbitalsExtrapolationOrder3<OrbitalsType>::extrapolate_orbitals(
    OrbitalsType** orbitals, OrbitalsType* new_orbitals)
{
    Control& ct = *(Control::instance());

    new_orbitals->assign(**orbitals);

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

    // do the extrapolation if previous orbitals exist (not at first step)

    if (orbitals_minus1_ != nullptr)
    {
        OrbitalsType tmp_orbitals_minus1("minus1", *new_orbitals, false);

        if (ct.verbose > 1 && onpe0)
            (*MPIdata::sout) << "Extrapolate orbitals using 3rd order scheme..."
                             << std::endl;
            // align orbitals_minus1 with new_orbitals
#ifdef MGMOL_USE_SCALAPACK
        if (use_dense_proj_mat)
        {
            dist_matrix::DistMatrix<DISTMATDTYPE> matQ("Q", ct.numst, ct.numst);
            dist_matrix::DistMatrix<DISTMATDTYPE> yyt(
                "yyt", ct.numst, ct.numst);

            // alignement
            orbitals_minus1_->computeGram(*new_orbitals, matQ);
            getProcrustesTransform(matQ, yyt);
            orbitals_minus1_->multiply_by_matrix(matQ);

            // compute delta Phi
            tmp_orbitals_minus1.assign(*orbitals_minus1_);
            tmp_orbitals_minus1.axpy((ORBDTYPE)-1., *new_orbitals);
            tmp_orbitals_minus1.multiply_by_matrix(yyt);

            if (orbitals_minus2_ != nullptr)
            {
                // alignement
                orbitals_minus2_->computeGram(*new_orbitals, matQ);
                getProcrustesTransform(matQ, yyt);
                orbitals_minus2_->multiply_by_matrix(matQ);

                // compute delta Phi
                orbitals_minus2_->axpy((ORBDTYPE)-1., *orbitals_minus1_);
                orbitals_minus2_->multiply_by_matrix(yyt);
            }
        }
        else
#endif
        {
            tmp_orbitals_minus1.assign(*orbitals_minus1_);
            if (orbitals_minus2_ != nullptr)
                orbitals_minus2_->axpy((ORBDTYPE)-1., *orbitals_minus1_);
            if (ct.verbose > 1 && onpe0)
                (*MPIdata::sout)
                    << "Compute tmp_orbitals_minus1..." << std::endl;
            tmp_orbitals_minus1.axpy((ORBDTYPE)-1., *new_orbitals);
        }

        if (orbitals_minus2_ != nullptr)
        {
            new_orbitals->axpy((ORBDTYPE)-2., tmp_orbitals_minus1);
            new_orbitals->axpy((ORBDTYPE)1., *orbitals_minus2_);

            delete orbitals_minus2_;
        }
        else
        {
            if (ct.verbose > 1 && onpe0)
                (*MPIdata::sout) << "Extrapolate orbitals using 2nd order "
                                    "scheme only for this step..."
                                 << std::endl;
            new_orbitals->axpy((ORBDTYPE)-1., tmp_orbitals_minus1);
        }

        orbitals_minus2_ = orbitals_minus1_;
    }

    orbitals_minus1_ = *orbitals;

    *orbitals = new_orbitals;

    (*orbitals)->incrementIterativeIndex();

    if (ct.isLocMode())
    {
        (*orbitals)->normalize();
        (*orbitals)->applyMask();
    }
    else
    {
        (*orbitals)->orthonormalizeLoewdin();
    }
}

template class OrbitalsExtrapolationOrder3<LocGridOrbitals<ORBDTYPE>>;
template class OrbitalsExtrapolationOrder3<ExtendedGridOrbitals<ORBDTYPE>>;
