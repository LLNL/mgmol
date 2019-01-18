// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: SubMatrices.h 7 2011-03-08 18:50:31Z jeanluc $
// Class SubMatrices
// container for m_ matrices n_*n_ corresponding to submatrices contributions
#ifndef SubMatrices_H
#define SubMatrices_H

#if USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;

#include "DistMatrix.h"
#include "SubMatricesIndexing.h"
#include "Timer.h"
// template <class T> class SubMatricesIndexing;

namespace dist_matrix
{

template <class T>
class SubMatrices
{
private:
    static Timer gather_tm_;
    static Timer gather_comp_tm_;
    static Timer gather_comm_tm_;

    string object_name_;

    // MPI
    const MPI_Comm comm_;
    int npes_;
    int mype_;
    int npes_distmat_;

    const SubMatricesIndexing<T>& submat_indexing_;

    // "local" data
    int n_; // size of local matrices
    int nb_local_matrices_; // number of matrices on this PE
    vector<T*> local_data_; // "local" data
    T* storage_;

    bool send_;

public:
    static Timer gather_tm() { return gather_tm_; }

    static Timer gather_comp_tm() { return gather_comp_tm_; }

    static Timer gather_comm_tm() { return gather_comm_tm_; }

    SubMatrices<T>(const string& name, const vector<vector<int>>& indexes,
        MPI_Comm comm, const DistMatrix<T>& mat, const SubMatricesIndexing<T>&);

    ~SubMatrices<T>() { delete[] storage_; }

    void scal(const double alpha);
    void gather(const DistMatrix<T>& mat);
    void print(ostream&);
    void get_array(const int im, vector<T>& val);

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
