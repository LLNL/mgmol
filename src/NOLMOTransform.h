// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

////////////////////////////////////////////////////////////////////////////////
//
// NOLMOTransform.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id$

#ifndef NOLMOTRANSFORM_H
#define NOLMOTRANSFORM_H

#include "OrbitalsTransform.h"
#include "Vector3D.h"

#include <string>
#include <vector>

class NOLMOTransform : public OrbitalsTransform
{
private:
    // cosine and sine matrices for x,y,z directions
    std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>*> b_;

    bool set_positions_;
    std::vector<std::vector<double>> ara0_;

    // lower cholesky decomposition transposed
    std::vector<DISTMATDTYPE> lst_;

    double nolmo(int maxsweep, double tol);
    double nolmo_fixedCenters(const int maxsweep, const double tol);

    double getFunctionalValue(const std::vector<std::vector<double>>& ara,
        const std::vector<std::vector<double>>& aba,
        const std::vector<double>& inv_alphan);
    double getFunctionalValueFixedCenters(
        const std::vector<std::vector<double>>& ara,
        const std::vector<std::vector<double>>& aba,
        const std::vector<double>& inv_alphan);
    void getNorms(std::vector<double>&, std::vector<double>&);
    double getDt();

public:
    std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>*>& b(void) { return b_; }

    dist_matrix::DistMatrix<DISTMATDTYPE>& b(const int i)
    {
        assert(i < 2 * NDIM && i >= 0);
        return *b_[i];
    }

    double spread2(int i, int j) const override;

    void init_transform(const dist_matrix::DistMatrix<DISTMATDTYPE>& ls,
        const bool reset_positions = false);

    // compute NOLMO transform
    void compute_transform(const int maxsweep, const double tol) override;

    NOLMOTransform(const int nst, const Vector3D& origin, const Vector3D& ll)
        : OrbitalsTransform(nst, origin, ll)
    {
        a_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(
            "matrixA", *bcr_, nst_, nst_, nst_, bsize_);
        const int n = a_->nloc(); // nb columns to compute on this PE
        if (bcr_->mycol() > -1)
            for (int i = 0; i < n; i++)
            {
                a_->setval(i * nst_ + i + offset_, 1.);
            }

        std::vector<std::string> names(2 * NDIM);
        names[0] = "matrixR_0";
        names[1] = "matrixR_1";
        names[2] = "matrixR_2";
        names[3] = "matrixR_3";
#if NDIM > 2
        names[4] = "matrixR_4";
        names[5] = "matrixR_5";
#endif
        for (int k = 0; k < 2 * NDIM; k++)
        {
            r_[k] = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                names[k], *bcr_, nst_, nst_, nst_, bsize_);
        }

        b_.resize(2 * NDIM);
        ara0_.resize(2 * NDIM);
        int n2   = nst_ * nst_;
        names[0] = "matrixB_0";
        names[1] = "matrixB_1";
        names[2] = "matrixB_2";
        names[3] = "matrixB_3";
#if NDIM > 2
        names[4] = "matrixB_4";
        names[5] = "matrixB_5";
#endif
        for (int k = 0; k < 2 * NDIM; k++)
        {
            b_[k] = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                names[k], *bcr_, nst_, nst_, nst_, bsize_);
        }
        for (int k = 0; k < 2 * NDIM; k++)
        {
            ara0_[k].resize(lnst_);
        }
        lst_.assign(n2, 0.);
        set_positions_ = false;
    };

    ~NOLMOTransform() override
    {
        delete a_;
        for (int k = 0; k < 2 * NDIM; k++)
        {
            assert(b_[k] != NULL);
            delete b_[k];
        }
    }

    void gatherTransformMat();
};
#endif
