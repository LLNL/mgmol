// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SquareLocalMatrices.h"

using namespace std;

template <class T>
SquareLocalMatrices<T>::SquareLocalMatrices(const int subdiv, const int m)
    : LocalMatrices<T>(subdiv, m, m){};

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

//set matrix elements to zero in rows/columns
//not associated with any orbital
template <class T>
void SquareLocalMatrices<T>::applySymmetricMask(const vector< vector<int>>& gids)
{
    const int m = LocalMatrices<T>::m_;
    const int n = LocalMatrices<T>::n_;

    for (short iloc = 0; iloc < LocalMatrices<T>::subdiv_; iloc++)
    {
        T* mat = LocalMatrices<T>::getSubMatrix(iloc);
        const vector<int>& loc_gids(gids[iloc]);

        for (int j = 0; j < n; j++)
        {
            const bool jvalid = (loc_gids[j] != -1);
            const int offset = j*m;

            for (int i = 0; i < m; i++)
            {
                const bool ivalid = (loc_gids[i] != -1);
                if ( !(ivalid && jvalid) )mat[offset + i] = 0.;
            }
        }
    }
}

template class SquareLocalMatrices<double>;
template class SquareLocalMatrices<float>;
