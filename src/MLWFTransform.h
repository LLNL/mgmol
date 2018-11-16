// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

////////////////////////////////////////////////////////////////////////////////
//
// MLWFTransform.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id$
// Adapted from  MLWFTransform.h,v 1.4 2002/06/28 18:05:58 fgygi Exp

#ifndef MLWFTRANSFORM_H
#define MLWFTRANSFORM_H

#include "OrbitalsTransform.h"
#include "Vector3D.h"

//#if 0
class MLWFTransform : public OrbitalsTransform
{
private:
    // total number of columns used in matrix, including dummy
    int nstcol_;

    double jade(int maxsweep, double tol);

public:
    double spread2(int i, int j) const;

    void setia(std::vector<int>&);

    // compute MLWF transform
    void compute_transform(const int maxsweep, const double tol);

    MLWFTransform(const int nst, const Vector3D& origin, const Vector3D& ll)
        : OrbitalsTransform(nst, origin, ll)
    {
        const int npcol = bcr_->npcol();
        nstcol_         = bsize_ * npcol;

        if (nst > 0)
        {
            a_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                "matrixA", *bcr_, nst_, nstcol_, nst_, bsize_);
            const int n = a_->nloc(); // nb columns to compute on this PE
            if (bcr_->mycol() > -1)
                for (int i = 0; i < n; i++)
                {
                    if ((offset_ + i) < nst_)
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
                    names[k], *bcr_, nst_, nstcol_, nst_, bsize_);
            }
        }
    };

    ~MLWFTransform() {}

    void printTransform();
};
//#endif
#endif
