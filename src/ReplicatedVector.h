// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef MGMOL_REPLICATEDVECTOR_H
#define MGMOL_REPLICATEDVECTOR_H

#include "ReplicatedMatrix.h"

#include <vector>

class ReplicatedVector
{
    int dim_;

    std::unique_ptr<double, void (*)(double*)> data_;

public:
    ReplicatedVector(const int n);
    ReplicatedVector(const ReplicatedVector&);
    ReplicatedVector(const std::vector<double>&);
    ReplicatedVector& operator=(const ReplicatedVector&);
    double* data() { return data_.get(); }
    void clear();
    double dot(const ReplicatedVector& v);
    double nrm2();
    void copyDataToVector(std::vector<double>& dst) const;
    void assign(const std::vector<double>& src);
    void axpy(const double alpha, const ReplicatedVector&);
    void gemv(const char trans, const double alpha, const ReplicatedMatrix&,
        const ReplicatedVector&, const double beta);
    void gemm(const char transa, const char transb, const double alpha,
        const ReplicatedMatrix&, const ReplicatedVector&, const double beta);
};

#endif
