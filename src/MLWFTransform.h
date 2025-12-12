// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_MLWFTRANSFORM_H
#define MGMOL_MLWFTRANSFORM_H

#include "OrbitalsTransform.h"
#include "Vector3D.h"

class MLWFTransform : public OrbitalsTransform
{
private:
    // total number of columns used in matrix, including dummy
    int nstcol_;

public:
    double spread2(int i, int j) const override;

    void setia(std::vector<int>&);

    // compute MLWF transform
    void compute_transform(const int maxsweep, const double tol) override;

    MLWFTransform(const int nst, const Vector3D& origin, const Vector3D& ll);

    ~MLWFTransform() override { }

    void printTransform();
};
#endif
