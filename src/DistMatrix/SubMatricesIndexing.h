// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_SubMatricesIndexing_H
#define MGMOL_SubMatricesIndexing_H

#include "DistMatrix.h"

#include <mpi.h>

#include <map>
#include <vector>

namespace dist_matrix
{

template <class T>
class SubMatricesIndexing
{
private:
    // MPI
    const MPI_Comm comm_;
    MPI_Comm distmat_active_comm_;
    int npes_;
    int mype_;
    int npes_distmat_;

    // arrays for MPI data transfer
    std::vector<int> double_local_indexes_;
    std::vector<int> my_sizes_;
    std::vector<int> remote_double_indexes_;
    std::vector<int> remote_sizes_;

    // maping to know location of index in array
    std::vector<std::map<int, short>>
        map_indexes_; // map function index <-> color
    std::vector<std::vector<int>> vector_indexes_;

    int nb_local_matrices_; // number of matrices on this PE
    const std::vector<std::vector<int>>
        local_indexes_; // indexes of "local" data

public:
    SubMatricesIndexing<T>(const std::vector<std::vector<int>>& indexes,
        MPI_Comm comm, const DistMatrix<T>& mat);

    int getLocalIndex(const int imat, const int index) const
    {
        return local_indexes_[imat][index];
    }
    int getDoubleIndex(const int i) const { return double_local_indexes_[i]; }

    int getRemoteSize(const int i) const
    {
        assert(i < (int)remote_sizes_.size());
        return remote_sizes_[i];
    }
    int getMySize(const int pe_distmat) const
    {
        if (pe_distmat < npes_distmat_)
            return my_sizes_[pe_distmat];
        else
            return 0;
    }
    int getRemoteDoubleIndex(const int i) const
    {
        assert(i < static_cast<int>(remote_double_indexes_.size()));
        return remote_double_indexes_[i];
    }
    int getVectorIndex(const int imat, const int index) const
    {
        return vector_indexes_[imat][index];
    }
    const std::map<int, short>& mapIndexes(const int imat) const
    {
        return map_indexes_[imat];
    }
};
}
#endif
