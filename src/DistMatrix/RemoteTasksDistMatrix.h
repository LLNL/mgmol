// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_REMOTETASKSDISTMATRIX_H
#define MGMOL_REMOTETASKSDISTMATRIX_H

#include "DistMatrix.h"

#include <cassert>
#include <climits>

#include <mpi.h>

namespace dist_matrix
{

template <class T>
class RemoteTasksDistMatrix
{
private:
    short my_task_index_;
    std::vector<short> remote_tasks_indexes_;

public:
    RemoteTasksDistMatrix(DistMatrix<T>& mat)
    {
        int npes            = 1;
        const MPI_Comm comm = mat.comm_global();
        MPI_Comm_size(comm, &npes);
        const short mypr  = mat.myrow();
        const short mypc  = mat.mycol();
        const short nprow = mat.nprow();
        assert(mypr + mypc * nprow < SHRT_MAX);
        my_task_index_ = mypr + mypc * nprow;

        remote_tasks_indexes_.resize(npes);
        if (npes > 1)
        {
            MPI_Allgather(&my_task_index_, 1, MPI_SHORT,
                &remote_tasks_indexes_[0], 1, MPI_SHORT, comm);
        }
        else
        {
            remote_tasks_indexes_[0] = my_task_index_;
        }
    }

    short getRemoteTask(const int pe) const
    {
        return remote_tasks_indexes_[pe];
    }

    bool isRemoteTaskActive(const int pe) const
    {
        return (remote_tasks_indexes_[pe] > -1);
    }

    short getMyTaskIndex() const { return my_task_index_; }

    bool isMyTaskActive() const { return (my_task_index_ > -1); }
};
}

#endif
