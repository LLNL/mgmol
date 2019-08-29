// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PROJECTED_MATRICES_MEHR_H
#define PROJECTED_MATRICES_MEHR_H

#include "DistMatrix.h"
#include "ProjectedMatrices.h"

#include <iostream>

class ProjectedMatricesMehrstellen : public ProjectedMatrices
{
    dist_matrix::DistMatrix<DISTMATDTYPE>* matB_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* invB_;

    void printB(std::ostream& os) const
    {
        assert(matB_ != 0);
        if (onpe0) os << " Matrix B" << std::endl;
        matB_->print(os, 0, 0, NPRINT_ROWS_AND_COLS, NPRINT_ROWS_AND_COLS);
    }

    void printInvB(std::ostream& os) const
    {
        assert(invB_ != 0);

        if (onpe0) os << " Matrix invB" << std::endl;
        invB_->print(os, 0, 0, NPRINT_ROWS_AND_COLS, NPRINT_ROWS_AND_COLS);
    }
    ProjectedMatricesMehrstellen(const ProjectedMatricesMehrstellen& pm);

public:
    ProjectedMatricesMehrstellen(const int ndim, const bool with_spin);
    ~ProjectedMatricesMehrstellen() override;

    void computeInvB() override;
    void rotateAll(const dist_matrix::DistMatrix<DISTMATDTYPE>& rotation_matrix,
        const bool flag_eigen) override;

    void initializeMatB(const SquareLocalMatrices<MATDTYPE>& ss) override
    {
        LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();

        sl2dm->accumulate(ss, *work_, dim_);

        *matB_ = *work_;
    }

    void updateTheta() override
    {
        // theta = invB * Hij
        theta_->symm('l', 'l', 1., *invB_, *matH_, 0.);
    }

    void updateHB() override
    {
        // if( onpe0 )
        //    (*MPIdata::sout)<<"ProjectedMatrices::updateHB()..."<<endl;
        matHB_->symm('l', 'l', 1., gm_->getMatrix(), *theta_, 0.);
        dist_matrix::DistMatrix<DISTMATDTYPE> work_dis(*matHB_);
        matHB_->transpose(0.5, work_dis, 0.5);
    }

    void printMatrices(std::ostream& os) const override
    {
        printS(os);
        printB(os);
        printInvB(os);
        printH(os);
        printTheta(os);
        printHB(os);
    }
};

#endif
