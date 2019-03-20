// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ProjectedMatricesMehrstellen.h"

#include "tools.h"

using namespace std;

ProjectedMatricesMehrstellen::ProjectedMatricesMehrstellen(
    const int ndim, const bool with_spin)
    : ProjectedMatrices(ndim, with_spin)
{
    if (dim_ > 0)
    {
        matB_ = new dist_matrix::DistMatrix<DISTMATDTYPE>("B", ndim, ndim);
        invB_ = new dist_matrix::DistMatrix<DISTMATDTYPE>("invB", ndim, ndim);
        matB_->identity();
        invB_->identity();
    }
}

ProjectedMatricesMehrstellen::~ProjectedMatricesMehrstellen()
{
    if (dim_ > 0)
    {
        delete matB_;
        matB_ = 0;
        delete invB_;
        invB_ = 0;
    }
}

void ProjectedMatricesMehrstellen::computeInvB()
{
    assert(matB_ != 0);
    assert(invB_ != 0);
    {

#ifdef DEBUG
        if (onpe0)
            (*MPIdata::sout) << "ProjectedMatrices::computeInvB()" << endl;
#endif
        *invB_    = *matB_;
        int info1 = invB_->potrf('l');
#ifdef USE_MPI
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.bcast(&info1, 1);
#endif
        if (info1 != 0)
        {
            if (onpe0) (*MPIdata::serr) << "Matrix: " << matB_->name() << endl;
            matB_->printMM((*MPIdata::serr));
        }
        int info2 = invB_->potri('l');
#ifdef USE_MPI
        mmpi.bcast(&info2, 1);
#endif
        if (info2 != 0 && onpe0)
        {
            (*MPIdata::serr) << "Matrix: " << matB_->name() << endl;
            matB_->printMM((*MPIdata::serr));
        }
    }
}

void ProjectedMatricesMehrstellen::rotateAll(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& rotation_matrix,
    const bool flag_eigen)
{
    // S -> U^T S U
    // rotate overlap and l_s
    if (flag_eigen)
    {
        gm_->set2Id(-1);
    }
    else
    {
        gm_->rotateAll(rotation_matrix);
    }
    //(*MPIdata::sout)<<"matS"<<endl;
    // matS_->print((*MPIdata::sout),0,0,5,5);

    // rotate matH_
    rotateSym(*matH_, rotation_matrix, *work_);

    // rotate matB_
    rotateSym(*matB_, rotation_matrix, *work_);

    computeInvB();

    // theta = invB * matH_
    updateTheta();

    updateHB();

    dm_->rotate(rotation_matrix, flag_eigen);
}
