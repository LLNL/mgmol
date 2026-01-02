// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ProjectedMatricesMehrstellen.h"
#include "ReplicatedMatrix.h"

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrix.h"
#include "DistMatrixTools.h"
#endif

template <class MatrixType>
ProjectedMatricesMehrstellen<MatrixType>::ProjectedMatricesMehrstellen(
    const int ndim, const bool with_spin, const double width)
    : ProjectedMatrices<MatrixType>(ndim, with_spin, width)
{
    assert(ndim > 0);

    matB_ = new MatrixType("B", ndim, ndim);
    invB_ = new MatrixType("invB", ndim, ndim);
    matB_->identity();
    invB_->identity();
}

template <class MatrixType>
ProjectedMatricesMehrstellen<MatrixType>::~ProjectedMatricesMehrstellen()
{
    delete matB_;
    matB_ = nullptr;
    delete invB_;
    invB_ = nullptr;
}

template <class MatrixType>
void ProjectedMatricesMehrstellen<MatrixType>::computeInvB()
{
    assert(matB_ != nullptr);
    assert(invB_ != nullptr);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
#ifdef DEBUG
    if (mmpi.instancePE0())
        (*MPIdata::sout) << "ProjectedMatrices::computeInvB()" << std::endl;
#endif
    *invB_    = *matB_;
    int info1 = invB_->potrf('l');
    mmpi.bcast(&info1, 1);
    if (info1 != 0)
    {
        if (mmpi.instancePE0())
            (*MPIdata::serr) << "Matrix: " << matB_->name() << std::endl;
        matB_->printMM((*MPIdata::serr));
    }
    int info2 = invB_->potri('l');
    mmpi.bcast(&info2, 1);
    if (info2 != 0 && mmpi.instancePE0())
    {
        (*MPIdata::serr) << "Matrix: " << matB_->name() << std::endl;
        matB_->printMM((*MPIdata::serr));
    }
}

template <class MatrixType>
void ProjectedMatricesMehrstellen<MatrixType>::rotateAll(
    const MatrixType& rotation_matrix, const bool flag_eigen)
{
    // S -> U^T S U
    // rotate overlap and l_s
    if (flag_eigen)
    {
        ProjectedMatrices<MatrixType>::gm_->set2Id(-1);
    }
    else
    {
        ProjectedMatrices<MatrixType>::gm_->rotateAll(rotation_matrix);
    }
    //(*MPIdata::sout)<<"matS"<<endl;
    // matS_->print((*MPIdata::sout),0,0,5,5);

    // rotate matH_
    rotateSym(*ProjectedMatrices<MatrixType>::matH_, rotation_matrix,
        *ProjectedMatrices<MatrixType>::work_);

    // rotate matB_
    rotateSym(*matB_, rotation_matrix, *ProjectedMatrices<MatrixType>::work_);

    computeInvB();

    // theta = invB * matH_
    updateTheta();

    updateHB();

    ProjectedMatrices<MatrixType>::dm_->rotate(rotation_matrix, flag_eigen);
}

#ifdef MGMOL_USE_SCALAPACK
template class ProjectedMatricesMehrstellen<
    dist_matrix::DistMatrix<DISTMATDTYPE>>;
#endif
template class ProjectedMatricesMehrstellen<ReplicatedMatrix>;
