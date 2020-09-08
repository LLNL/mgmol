// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SquareLocalMatrices.h"

template <class T>
SquareLocalMatrices<T>::SquareLocalMatrices(const int subdiv, const int m)
    : LocalMatrices<T>(subdiv, m, m)
{
}

template <class T>
void SquareLocalMatrices<T>::fillUpperWithLower()
{
    int m = LocalMatrices<T>::m_;

    for (short iloc = 0; iloc < LocalMatrices<T>::subdiv_; iloc++)
    {
        T* ssiloc = LocalMatrices<T>::getSubMatrix(iloc);

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

template <class T>
void SquareLocalMatrices<T>::setDiagonal2Zero()
{
    int m = LocalMatrices<T>::m_;

    for (short iloc = 0; iloc < LocalMatrices<T>::subdiv_; iloc++)
    {
        T* ssiloc = LocalMatrices<T>::getSubMatrix(iloc);
        for (int i = 0; i < m; i++)
        {
            ssiloc[i + m * i] = 0.;
        }
    }
}

template <class T>
void SquareLocalMatrices<T>::transpose()
{
    int m = LocalMatrices<T>::m_;

    for (short iloc = 0; iloc < LocalMatrices<T>::subdiv_; iloc++)
    {
        T* ssiloc = LocalMatrices<T>::getSubMatrix(iloc);
        MemorySpace::assert_is_host_ptr(ssiloc);

        for (int i = 0; i < m; i++)
        {
            const int istart = m * i;
            for (int j = 0; j < i; j++)
            {
                T tmp              = ssiloc[i + m * j];
                ssiloc[i + m * j]  = ssiloc[istart + j];
                ssiloc[istart + j] = tmp;
            }
        }
    }
}

template <class T>
double SquareLocalMatrices<T>::computePartialTrace(
    const std::vector<int>& ids, const int iloc)
{
    assert(!ids.empty());

    int m = LocalMatrices<T>::m_;

    T* ssiloc = LocalMatrices<T>::getSubMatrix(iloc);
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

template <class T>
void SquareLocalMatrices<T>::shift(const T shift)
{
    const int m = LocalMatrices<T>::m_;

    for (short iloc = 0; iloc < LocalMatrices<T>::subdiv_; iloc++)
    {
        T* mat = LocalMatrices<T>::getSubMatrix(iloc);
        MemorySpace::assert_is_host_ptr(mat);

        for (int i = 0; i < m; i++)
        {
            mat[i + m * i] += shift;
        }
    }
}

// set matrix elements to zero in rows/columns
// not associated with any orbital
template <class T>
void SquareLocalMatrices<T>::applySymmetricMask(
    const std::vector<std::vector<int>>& gids)
{
    const int m = LocalMatrices<T>::m_;
    const int n = LocalMatrices<T>::n_;

    for (short iloc = 0; iloc < LocalMatrices<T>::subdiv_; iloc++)
    {
        T* mat = LocalMatrices<T>::getSubMatrix(iloc);
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
