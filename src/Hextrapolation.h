// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_HEXTRAPLOATION_H
#define MGMOL_HEXTRAPLOATION_H

#include "DistMatrix.h"
#include "MPIdata.h"
#include "tools.h"

#include <iostream>

class Hextrapolation
{
    int ndim_;

    dist_matrix::DistMatrix<DISTMATDTYPE>* h_minus1_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* h_minus2_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* new_h_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* tmp_h_minus1_;

    dist_matrix::DistMatrix<DISTMATDTYPE>* work_;

public:
    Hextrapolation(const int ndim) : ndim_(ndim)
    {
        work_ = new dist_matrix::DistMatrix<DISTMATDTYPE>("work", ndim, ndim);
    }

    ~Hextrapolation() { delete work_; }

    void initExtrapolationH(const dist_matrix::DistMatrix<DISTMATDTYPE>& matHB)
    {
        *new_h_ = matHB;
    }

    void extrapolateHorder2(const dist_matrix::DistMatrix<DISTMATDTYPE>& matQ,
        dist_matrix::DistMatrix<DISTMATDTYPE>& yyt, std::ostream& os)
    {
        rotateSym(*h_minus1_, matQ, *work_);

        h_minus1_->axpy(-1., *new_h_);
        work_->gemm('n', 'n', 1., *h_minus1_, yyt, 0.);
        h_minus1_->gemm('t', 'n', 1., yyt, *work_, 0.);
        if (onpe0) os << "delta H..." << std::endl;
        h_minus1_->print(os, 0, 0, 5, 5);
        new_h_->axpy(-1., *h_minus1_);
    }

    void updateHminus1tmp(const dist_matrix::DistMatrix<DISTMATDTYPE>& matQ,
        dist_matrix::DistMatrix<DISTMATDTYPE>& yyt, std::ostream& os)
    {
        if (h_minus1_ != nullptr)
        {
            rotateSym(*h_minus1_, matQ, *work_);
            *tmp_h_minus1_ = *h_minus1_;
            tmp_h_minus1_->axpy(-1., *new_h_);
            work_->gemm('n', 'n', 1., *tmp_h_minus1_, yyt, 0.);
            tmp_h_minus1_->gemm('t', 'n', 1., yyt, *work_, 0.);
            if (onpe0) os << "delta H..." << std::endl;
            tmp_h_minus1_->print(os, 0, 0, 5, 5);
        }
    }

    void updateHminus2(const dist_matrix::DistMatrix<DISTMATDTYPE>& matQ,
        dist_matrix::DistMatrix<DISTMATDTYPE>& yyt, std::ostream& os)
    {
        if (h_minus2_ != nullptr)
        {
            rotateSym(*h_minus2_, matQ, *work_);
            h_minus2_->axpy(-1., *h_minus1_);
            work_->gemm('n', 'n', 1., *h_minus2_, yyt, 0.);
            h_minus2_->gemm('t', 'n', 1., yyt, *work_, 0.);
            if (onpe0) os << "delta H..." << std::endl;
            h_minus2_->print(os, 0, 0, 5, 5);
        }
    }

    void extrapolateHorder3()
    {
        if (h_minus2_ != nullptr)
        {
            new_h_->axpy(-2., *tmp_h_minus1_);
            new_h_->axpy(1., *h_minus2_);
        }
        else
        {
            if (h_minus1_ != nullptr)
            {
                new_h_->axpy(-1., *tmp_h_minus1_);
            }
        }
    }

    void saveH(dist_matrix::DistMatrix<DISTMATDTYPE>& matHB)
    {
        matHB = *new_h_;
    }

    void updateHminus1(const dist_matrix::DistMatrix<DISTMATDTYPE>& matHB)
    {
        if (h_minus1_ == nullptr)
        {
            h_minus1_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                "H_minus1", ndim_, ndim_);
        }
        *h_minus1_ = matHB;
    }

    void updateHminus2(const dist_matrix::DistMatrix<DISTMATDTYPE>& matHB)
    {
        if (h_minus1_ == nullptr)
        {
            h_minus1_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                "H_minus1", ndim_, ndim_);
        }
        else
        {
            if (h_minus2_ == nullptr)
                h_minus2_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                    "H_minus2", ndim_, ndim_);
            *h_minus2_ = *h_minus1_;
        }
        *h_minus1_ = matHB;
    }
};

#endif
