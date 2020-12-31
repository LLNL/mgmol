// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SquareLocalMatrices.h"
#ifdef HAVE_MAGMA
#include "ReplicatedMatrix.h"
#endif

template <typename DataType, typename MemorySpaceType>
SquareLocalMatrices<DataType, MemorySpaceType>::SquareLocalMatrices(
    const int nmat, const int m)
    : LocalMatrices<DataType, MemorySpaceType>(nmat, m, m)
{
}

template <typename DataType, typename MemorySpaceType>
void SquareLocalMatrices<DataType, MemorySpaceType>::fillUpperWithLower()
{
    int m = LocalMatrices<DataType, MemorySpaceType>::m_;

    for (short iloc = 0; iloc < LocalMatrices<DataType, MemorySpaceType>::nmat_;
         iloc++)
    {
        DataType* ssiloc
            = LocalMatrices<DataType, MemorySpaceType>::getRawPtr(iloc);

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

template <typename DataType, typename MemorySpaceType>
void SquareLocalMatrices<DataType, MemorySpaceType>::setDiagonal2Zero()
{
    int m = LocalMatrices<DataType, MemorySpaceType>::m_;

    for (short iloc = 0; iloc < LocalMatrices<DataType, MemorySpaceType>::nmat_;
         iloc++)
    {
        DataType* ssiloc
            = LocalMatrices<DataType, MemorySpaceType>::getRawPtr(iloc);
        for (int i = 0; i < m; i++)
        {
            ssiloc[i + m * i] = 0.;
        }
    }
}

template <typename DataType, typename MemorySpaceType>
void SquareLocalMatrices<DataType, MemorySpaceType>::transpose()
{
    int m = LocalMatrices<DataType, MemorySpaceType>::m_;

    for (short iloc = 0; iloc < LocalMatrices<DataType, MemorySpaceType>::nmat_;
         iloc++)
    {
        DataType* ssiloc
            = LocalMatrices<DataType, MemorySpaceType>::getRawPtr(iloc);
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

template <typename DataType, typename MemorySpaceType>
double SquareLocalMatrices<DataType, MemorySpaceType>::computePartialTrace(
    const std::vector<int>& ids, const int iloc)
{
    assert(!ids.empty());

    int m = LocalMatrices<DataType, MemorySpaceType>::m_;

    DataType* ssiloc
        = LocalMatrices<DataType, MemorySpaceType>::getRawPtr(iloc);
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

template <typename DataType, typename MemorySpaceType>
void SquareLocalMatrices<DataType, MemorySpaceType>::shift(const DataType shift)
{
    const int m = LocalMatrices<DataType, MemorySpaceType>::m_;

    for (short iloc = 0; iloc < LocalMatrices<DataType, MemorySpaceType>::nmat_;
         iloc++)
    {
        DataType* mat
            = LocalMatrices<DataType, MemorySpaceType>::getRawPtr(iloc);
        MemorySpace::assert_is_host_ptr(mat);

        for (int i = 0; i < m; i++)
        {
            mat[i + m * i] += shift;
        }
    }
}

// set matrix elements to zero in rows/columns
// not associated with any orbital
template <typename DataType, typename MemorySpaceType>
void SquareLocalMatrices<DataType, MemorySpaceType>::applySymmetricMask(
    const std::vector<std::vector<int>>& gids)
{
    const int m = LocalMatrices<DataType, MemorySpaceType>::m_;
    const int n = LocalMatrices<DataType, MemorySpaceType>::n_;

    for (short iloc = 0; iloc < LocalMatrices<DataType, MemorySpaceType>::nmat_;
         iloc++)
    {
        DataType* mat
            = LocalMatrices<DataType, MemorySpaceType>::getRawPtr(iloc);
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

#ifdef HAVE_MAGMA
template <>
void SquareLocalMatrices<double, MemorySpace::Device>::assign(
    const ReplicatedMatrix& src, const int iloc)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    double* dst = LocalMatrices<double, MemorySpace::Device>::getRawPtr(iloc);
    magma_dcopymatrix(src.m(), src.m(), src.data(), src.ld(), dst, m_,
        magma_singleton.queue_);
}

template <>
void SquareLocalMatrices<double, MemorySpace::Device>::assign(
    const SquareLocalMatrices<double, MemorySpace::Host>& src)
{
    auto& magma_singleton = MagmaSingleton::get_magma_singleton();

    for (short iloc = 0; iloc < nmat_; iloc++)
    {
        double* dst
            = LocalMatrices<double, MemorySpace::Device>::getRawPtr(iloc);
        magma_dsetmatrix(src.m(), src.n(), src.getSubMatrix(), src.n(), dst, m_,
            magma_singleton.queue_);
    }
}
#endif

template class SquareLocalMatrices<double, MemorySpace::Host>;
template class SquareLocalMatrices<float, MemorySpace::Host>;
#ifdef HAVE_MAGMA
template class SquareLocalMatrices<double, MemorySpace::Device>;
#endif
