// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// container for m_ matrices n_*n_ corresponding to submatrices contributions
#ifndef MGMOL_SubMatrices_H
#define MGMOL_SubMatrices_H

#include <mpi.h>

#include "DistMatrix.h"
#include "SubMatricesIndexing.h"
#include "Timer.h"
// template <class T> class SubMatricesIndexing;

#include <cstring>
#include <iostream>
#include <vector>

typedef unsigned int type_displ;

namespace dist_matrix
{

template <class T>
class SubMatrices
{
private:
    static Timer gather_tm_;
    static Timer gather_comp_tm_;
    static Timer gather_comm_tm_;

    std::string object_name_;

    // MPI
    const MPI_Comm comm_;
    int npes_;
    int mype_;
    int npes_distmat_;

    const SubMatricesIndexing<T>& submat_indexing_;

    // "local" data
    int n_; // size of local matrices
    int nb_local_matrices_; // number of matrices on this PE
    std::vector<T*> local_data_; // "local" data
    T* storage_;

    bool send_;

    void sendrecv(T* sendbuf, T* recvbuf,
        const std::vector<type_displ>& my_displ,
        const std::vector<type_displ>& remote_displ);

    // templated MPI wrappers
    int internal_MPI_Irecv(
        T* buf, int count, int source, int tag, MPI_Request* request);
    int internal_MPI_Send(const T* buf, int count, int dest, int tag);

public:
    static void printTimers(std::ostream& os)
    {
        gather_tm_.print(os);
        gather_comp_tm_.print(os);
        gather_comm_tm_.print(os);
    }

    SubMatrices<T>(const std::string& name,
        const std::vector<std::vector<int>>& indexes, MPI_Comm comm,
        const DistMatrix<T>& mat, const SubMatricesIndexing<T>&);

    ~SubMatrices<T>() { delete[] storage_; }

    void scal(const double alpha);
    void gather(const DistMatrix<T>& mat);
    void print(std::ostream&) const;
    void get_array(const int im, std::vector<T>& val);

    T val(const int i, const int j, const int im) const
    {
        assert(im < nb_local_matrices_);
        assert(i >= 0);
        assert(j >= 0);
        assert(i < n_);
        assert(j < n_);
        return local_data_[im][i + n_ * j];
    }

    SubMatrices<T>& operator=(const SubMatrices<T>&);
};

template <class T>
Timer SubMatrices<T>::gather_tm_("SubMatrices<T>::gather");
template <class T>
Timer SubMatrices<T>::gather_comp_tm_("SubMatrices<T>::gather_comp");
template <class T>
Timer SubMatrices<T>::gather_comm_tm_("SubMatrices<T>::gather_comm");

} // namespace

#endif
