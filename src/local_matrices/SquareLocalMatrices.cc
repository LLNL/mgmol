// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SquareLocalMatrices.h"

template <class DataType>
SquareLocalMatrices<DataType>::SquareLocalMatrices(const int nmat, const int m)
    : LocalMatrices<DataType>(nmat, m, m)
{
}

template <class DataType>
void SquareLocalMatrices<DataType>::fillUpperWithLower()
{
    int m = LocalMatrices<DataType>::m_;

    for (short iloc = 0; iloc < LocalMatrices<DataType>::nmat_; iloc++)
    {
        DataType* ssiloc = LocalMatrices<DataType>::getRawPtr(iloc);

        for (int i = 0; i < m; i++)
        {
            const int istart = m * i;
            for (int j = 0; j < i; j++)
            {
                ssiloc[istart + j] = ssiloc[i + m * j];
            }
        }
    }
}

template <class DataType>
void SquareLocalMatrices<DataType>::setDiagonal2Zero()
{
    int m = LocalMatrices<DataType>::m_;

    for (short iloc = 0; iloc < LocalMatrices<DataType>::nmat_; iloc++)
    {
        DataType* ssiloc = LocalMatrices<DataType>::getRawPtr(iloc);
        for (int i = 0; i < m; i++)
        {
            ssiloc[i + m * i] = 0.;
        }
    }
}

template <class DataType>
void SquareLocalMatrices<DataType>::transpose()
{
    int m = LocalMatrices<DataType>::m_;

    for (short iloc = 0; iloc < LocalMatrices<DataType>::nmat_; iloc++)
    {
        DataType* ssiloc = LocalMatrices<DataType>::getRawPtr(iloc);
        MemorySpace::assert_is_host_ptr(ssiloc);

        for (int i = 0; i < m; i++)
        {
            const int istart = m * i;
            for (int j = 0; j < i; j++)
            {
                DataType tmp       = ssiloc[i + m * j];
                ssiloc[i + m * j]  = ssiloc[istart + j];
                ssiloc[istart + j] = tmp;
            }
        }
    }
}

template <class DataType>
double SquareLocalMatrices<DataType>::computePartialTrace(
    const std::vector<int>& ids, const int iloc)
{
    assert(!ids.empty());

    int m = LocalMatrices<DataType>::m_;

    DataType* ssiloc = LocalMatrices<DataType>::getRawPtr(iloc);
    MemorySpace::assert_is_host_ptr(ssiloc);

    double trace = 0.;
    const int n  = static_cast<int>(ids.size());
#pragma omp parallel for reduction(+ : trace)
    for (int i = 0; i < n; i++)
    {
        assert(ids[i] < m);
        trace += ssiloc[ids[i] * (m + 1)];
    }

    return trace;
}

template <class DataType>
void SquareLocalMatrices<DataType>::shift(const DataType shift)
{
    const int m = LocalMatrices<DataType>::m_;

    for (short iloc = 0; iloc < LocalMatrices<DataType>::nmat_; iloc++)
    {
        DataType* mat = LocalMatrices<DataType>::getRawPtr(iloc);
        MemorySpace::assert_is_host_ptr(mat);

        for (int i = 0; i < m; i++)
        {
            mat[i + m * i] += shift;
        }
    }
}

// set matrix elements to zero in rows/columns
// not associated with any orbital
template <class DataType>
void SquareLocalMatrices<DataType>::applySymmetricMask(
    const std::vector<std::vector<int>>& gids)
{
    const int m = LocalMatrices<DataType>::m_;
    const int n = LocalMatrices<DataType>::n_;

    for (short iloc = 0; iloc < LocalMatrices<DataType>::nmat_; iloc++)
    {
        DataType* mat = LocalMatrices<DataType>::getRawPtr(iloc);
        MemorySpace::assert_is_host_ptr(mat);
        const std::vector<int>& loc_gids(gids[iloc]);

        for (int j = 0; j < n; j++)
        {
            const bool jvalid = (loc_gids[j] != -1);
            const int offset  = j * m;

            for (int i = 0; i < m; i++)
            {
                const bool ivalid = (loc_gids[i] != -1);
                if (!(ivalid && jvalid)) mat[offset + i] = 0.;
            }
        }
    }
}

template class SquareLocalMatrices<double>;
template class SquareLocalMatrices<float>;
