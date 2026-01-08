// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DISTVECTOR_H
#define MGMOL_DISTVECTOR_H

#include "BlacsContext.h"
#include "DistMatrix.h"
#include "MGmol_MPI.h"

#include <string>
#include <vector>

namespace dist_matrix
{

template <class T>
class DistVector : public DistMatrix<T>
{
private:
    T operator[](int i) const { return DistMatrix<T>::val(i); }

public:
    DistVector(const std::string& name, const int m) : DistMatrix<T>(name, m, 1)
    {
    }

    DistVector(const int m) : DistMatrix<T>("noname", m, 1) { }

    DistVector(const std::vector<T>& v) : DistMatrix<T>("noname", v.size(), 1)
    {
        assign(v);
    }

    T* data() { return DistMatrix<T>::data(); }

    void assign(const std::vector<T>& v)
    {
        int m = (int)v.size();
        // Build vector fully located on PE0
        BlacsContext bc(DistMatrix<T>::comm_global(), 0, 1, 1);
        DistMatrix<T> tmp("tmp", bc, m, 1, m, 1);
        tmp.assign0(v);

        DistMatrix<T>::assign(tmp, 0, 0);
    }

    void swap(DistVector& v) { DistMatrix<T>::swap(v); }

    double nrm2() const
    {
        double dott = this->dot(*this);
        return sqrt(dott);
    }

    void normalize()
    {
        double norm2 = nrm2();
        DistMatrix<T>::scal(1. / norm2);
    }

    T scaledDiff2(const DistVector<T>& v, const double theta)
    {
        double diff2 = 0.;

        if (DistMatrix<T>::active() && v.nloc() > 0)
        {
            int m = v.mloc();
            for (int i = 0; i < m; i++)
            {
                double tmp = DistMatrix<T>::val_[i] - theta * v.val_[i];
                diff2 += tmp * tmp;
            }
        }
        double tmp      = diff2;
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.allreduce(&tmp, &diff2, 1, MPI_SUM);

        return diff2;
    }

    double dot(const DistVector<T>& v) const
    {
        double sum  = 0.;
        double tsum = 0.;
        if (DistMatrix<T>::active())
        {
            assert(v.size_ == DistMatrix<T>::size_);
            tsum = LinearAlgebraUtils<MemorySpace::Host>::MPdot(
                DistMatrix<T>::val_.size(), DistMatrix<T>::val_.data(),
                v.val_.data());
        }
#ifdef MGMOL_USE_SCALAPACK
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.allreduce(&tsum, &sum, 1, MPI_SUM);
#else
        sum = tsum;
#endif
        return sum;
    }
};
}
#endif
