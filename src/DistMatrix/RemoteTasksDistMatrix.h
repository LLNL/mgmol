// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: RemoteTasksDistMatrix.h 17 2012-08-15 17:11:26Z jeanluc $
#ifndef REMOTETASKSDISTMATRIX_H
#define REMOTETASKSDISTMATRIX_H

#include "DistMatrix.h"

#include <cassert>
#include <climits>

#if USE_MPI
#include <mpi.h>
#endif

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
        int npes = 1;
#if USE_MPI
        const MPI_Comm comm = mat.comm_global();
        MPI_Comm_size(comm, &npes);
#endif
        const short mypr  = mat.myrow();
        const short mypc  = mat.mycol();
        const short nprow = mat.nprow();
        assert(mypr + mypc * nprow < SHRT_MAX);
        my_task_index_ = mypr + mypc * nprow;

        remote_tasks_indexes_.resize(npes);
#if USE_MPI
        if (npes > 1)
        {
            MPI_Allgather(&my_task_index_, 1, MPI_SHORT,
                &remote_tasks_indexes_[0], 1, MPI_SHORT, comm);
        }
        else
#endif
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
