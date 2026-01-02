// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_PROJECTED_MATRICES_MEHR_H
#define MGMOL_PROJECTED_MATRICES_MEHR_H

#include "ProjectedMatrices.h"

#include <iostream>

template <class MatrixType>
class ProjectedMatricesMehrstellen : public ProjectedMatrices<MatrixType>
{
private:
    MatrixType* matB_;
    MatrixType* invB_;

    void printB(std::ostream& os) const
    {
        assert(matB_ != 0);
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        if (mmpi.instancePE0()) os << " Matrix B" << std::endl;
        matB_->print(os, 0, 0, NPRINT_ROWS_AND_COLS, NPRINT_ROWS_AND_COLS);
    }

    void printInvB(std::ostream& os) const
    {
        assert(invB_ != 0);
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        if (mmpi.instancePE0()) os << " Matrix invB" << std::endl;
        invB_->print(os, 0, 0, NPRINT_ROWS_AND_COLS, NPRINT_ROWS_AND_COLS);
    }
    ProjectedMatricesMehrstellen(const ProjectedMatricesMehrstellen& pm);

public:
    ProjectedMatricesMehrstellen(
        const int ndim, const bool with_spin, const double width);
    ~ProjectedMatricesMehrstellen() override;

    void computeInvB() override;
    void rotateAll(
        const MatrixType& rotation_matrix, const bool flag_eigen) override;

    void initializeMatB(
        const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& ss) override
    {
        ProjectedMatrices<MatrixType>::convert(ss, *matB_);
    }

    void updateTheta() override
    {
        // theta = invB * Hij
        ProjectedMatrices<MatrixType>::theta_->symm(
            'l', 'l', 1., *invB_, *ProjectedMatrices<MatrixType>::matH_, 0.);
    }

    void updateHB() override
    {
        ProjectedMatrices<MatrixType>::matHB_->symm('l', 'l', 1.,
            ProjectedMatrices<MatrixType>::gm_->getMatrix(),
            *ProjectedMatrices<MatrixType>::theta_, 0.);
        MatrixType work_dis(*ProjectedMatrices<MatrixType>::matHB_);
        ProjectedMatrices<MatrixType>::matHB_->transpose(0.5, work_dis, 0.5);
    }

    void printMatrices(std::ostream& os) const override
    {
        ProjectedMatrices<MatrixType>::printS(os);
        printB(os);
        printInvB(os);
        ProjectedMatrices<MatrixType>::printH(os);
        ProjectedMatrices<MatrixType>::printTheta(os);
        ProjectedMatrices<MatrixType>::printHB(os);
    }
};

#endif
