// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "OrbitalsExtrapolationOrder2.h"
#include "Control.h"
#include "DistMatrixTools.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "ProjectedMatrices.h"

#define EXTRAPOLATE_H 1
using namespace std;

template <class T>
void OrbitalsExtrapolationOrder2<T>::extrapolate_orbitals(
    T** orbitals, T* new_orbitals)
{
    Control& ct = *(Control::instance());

    bool use_dense_proj_mat = false;
#if EXTRAPOLATE_H
    ProjectedMatricesInterface* proj_matrices = (*orbitals)->getProjMatrices();
    ProjectedMatrices* projmat                = 0;
    if (ct.OuterSolver() != OuterSolverType::ABPG &&
        ct.OuterSolver() != OuterSolverType::NLCG )
    {
        projmat = dynamic_cast<ProjectedMatrices*>(proj_matrices);
        if( projmat )use_dense_proj_mat = true;
    }
#endif

    new_orbitals->assign(**orbitals);

    extrapolated_H_ = false;
    // do the extrapolation if previous orbitals exist (not at first step)
    {
        if (OrbitalsExtrapolation<T>::orbitals_minus1_ != 0)
        {
            T* orbitals_minus1 = OrbitalsExtrapolation<T>::orbitals_minus1_;
            if (ct.verbose > 1 && onpe0)
                (*MPIdata::sout) << "Extrapolate orbitals order 2..." << endl;

            // align orbitals_minus1_ with new_orbitals
            if (use_dense_proj_mat)
            {
                dist_matrix::DistMatrix<DISTMATDTYPE> matQ(
                    "Q", ct.numst, ct.numst);
#if 0
                T tmp(*orbitals_minus1_);
                tmp.axpy(-1.,*new_orbitals);
                tmp.computeGram(matQ);
                double normQ=matQ.trace();
                if( onpe0 )
                    (*MPIdata::sout)<<"||Phi_old-Phi_new|| before Procrustes = "<<normQ<<endl;
#endif
                orbitals_minus1->computeGram(*new_orbitals, matQ);

                dist_matrix::DistMatrix<DISTMATDTYPE> yyt(
                    "yyt", ct.numst, ct.numst);
                getProcrustesTransform(matQ, yyt);

                orbitals_minus1->multiply_by_matrix(matQ);
                orbitals_minus1->axpy(-1., *new_orbitals);
                orbitals_minus1->multiply_by_matrix(yyt);

#if EXTRAPOLATE_H
                if (ct.verbose > 2 && onpe0)
                    (*MPIdata::sout) << "Extrapolate H..." << endl;
                OrbitalsExtrapolation<T>::hextrapol_->extrapolateHorder2(
                    matQ, yyt, (*MPIdata::sout));
                extrapolated_H_ = true;
#endif

#if 0          
                tmp.assign(*orbitals_minus1_);
                //tmp.axpy(-1.,*new_orbitals);
                tmp.computeGram(matQ);
                normQ=matQ.trace();
                if( onpe0 )
                    (*MPIdata::sout)<<"||Phi_old-Phi_new|| after Procrustes = "<<normQ<<endl;
#endif
            }
            else
            { // !use_dense_proj_mat

                new_orbitals->scal(2.);
            }
            new_orbitals->axpy(-1., *orbitals_minus1);

            delete orbitals_minus1;
        }

        // save data for next extrapolation
        if (ct.verbose > 1 && onpe0)
            (*MPIdata::sout)
                << "Set orbitals_minus1_ to values of orbitals" << endl;
        OrbitalsExtrapolation<T>::orbitals_minus1_ = *orbitals;

        if (use_dense_proj_mat)
        {
#if EXTRAPOLATE_H
            hextrapol_->updateHminus1( projmat->getMatHB() );
#endif
        }
    }

    *orbitals = new_orbitals;
#if EXTRAPOLATE_H
    if (use_dense_proj_mat)
    {
        hextrapol_->saveH(projmat->getMatHB());
    }
#endif

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
        (*orbitals)->orthonormalizeLoewdin();
    }
}

template class OrbitalsExtrapolationOrder2<LocGridOrbitals>;
template class OrbitalsExtrapolationOrder2<ExtendedGridOrbitals>;
