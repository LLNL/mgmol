// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_KBPSIMATRIX_H
#define MGMOL_KBPSIMATRIX_H

#include <cassert>
#include <limits.h>
#include <set>
#include <vector>

#include "GridFunc.h"
#include "KBPsiMatrixInterface.h"
#include "Lap.h"
#include "LocGridOrbitals.h"
#include "Timer.h"

class Ions;
class Ion;
class ProjectedMatrices;
class DensityMatrix;

#if (UINT_MAX >= 4294967295)

typedef unsigned int INDEX_TYPE;
#define MPI_INDEX_TYPE MPI_UNSIGNED

#else

typedef unsigned long INDEX_TYPE;
#define MPI_INDEX_TYPE MPI_UNSIGNED_LONG

#endif

const int bit_shift         = 16;
const INDEX_TYPE modulo_int = (1 << bit_shift) - 1;

class KBPsiMatrix : public KBPsiMatrixInterface
{
    static Timer allreduce_tm_;
    static Timer global_sum_tm_;
    static Timer compute_kbpsi_tm_;
    static Timer computeHvnlMatrix_tm_;
    static Timer allGatherNonzeroElements_tm_;

    // storage arrays for KB states projections
    double* storage_;
    // ... and their size
    int size_kbpsi_;

    int numst_;
    int count_proj_subdomain_;

    // operator used to compute r.h.s
    // usually identity, unless Mehrstellen scheme is used
    const pb::Lap<ORBDTYPE>* lapop_;

    bool is_setup_;

    // pointers for columns of data corresponding to each ion
    // in array_kbpsi_
    vector<double*> ptr_kbpsi_;
    vector<double*> ptr_kbBpsi_;

    static set<INDEX_TYPE> nonzero_elements_;

    MPI_Comm comm_dir_[3];

    static dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>*
        remote_tasks_DistMatrix_;

    void addKBPsi(const int gid, const int st, const double val)
    {
        assert(gid < (int)ptr_kbpsi_.size());
        assert(st < numst_);
        assert(ptr_kbpsi_[gid] != 0);
        ptr_kbpsi_[gid][st] += val;
    }
    void addKBBPsi(const int gid, const int st, const double val)
    {
        assert(gid < (int)ptr_kbBpsi_.size());
        assert(st < numst_);
        assert(ptr_kbBpsi_[gid] != 0);
        ptr_kbBpsi_[gid][st] += val;
    }

    void allGatherNonzeroElements(MPI_Comm comm, const int ntasks);

    void clear();

    INDEX_TYPE index(const unsigned int st, const unsigned int ion) const
    {
        return (st << bit_shift) + ion + 1;
    }
    int index2st(const INDEX_TYPE index) const
    {
        return (int)(index >> bit_shift);
    }
    int index2ion(const INDEX_TYPE index) const
    {
        return (int)(index & modulo_int) - 1;
    }
    void reset()
    {
        if (storage_ != NULL)
            memset(storage_, 0, 2 * size_kbpsi_ * sizeof(double));
    }
    void computeSetNonZeroElements(
        const Ions& ions, const LocGridOrbitals& orbitals);
    void allocate(const Ions& ions, const int num_st);
    void globalSumKBpsi(const Ions& ions);

    void computeHvnlMatrix(const KBPsiMatrix* const kbpsi, const Ion&,
        dist_matrix::SparseDistMatrix<DISTMATDTYPE>&) const;
    void computeHvnlMatrix(const KBPsiMatrix* const kbpsi2, const Ion& ion,
        ProjectedMatricesInterface* proj_matrices) const;

    void getPsiKBPsiSym(
        const Ion& ion, dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sm);
    void getPsiKBPsiSym(
        const Ions& ions, dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sm);
    void getPsiKBPsiSym(
        const Ions& ions, dist_matrix::DistMatrix<DISTMATDTYPE>& Aij);

    void computeKBpsi(Ions& ions, LocGridOrbitals& orbitals,
        const int first_color, const int nb_colors, const bool flag);
    void computeKBpsi(
        Ions& ions, pb::GridFunc<ORBDTYPE>*, const int, const bool flag);

public:
    KBPsiMatrix(pb::Lap<ORBDTYPE>* lapop);

    ~KBPsiMatrix();

    void registerRemoteTasksDistMatrix(
        dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>*
            remote_tasks_DistMatrix)
    {
        assert(remote_tasks_DistMatrix != 0);
        remote_tasks_DistMatrix_ = remote_tasks_DistMatrix;
    }

    void printTimers(ostream& os);

    double getEvnl(const Ions& ions, LocGridOrbitals& orbitals,
        ProjectedMatrices* proj_matrices);
    double getValIonState(const int gid, const int st) const
    {
        assert(st < numst_);
        assert(gid < (int)ptr_kbpsi_.size());
        return ptr_kbpsi_[gid][st];
    }
    void scaleWithKBcoeff(const Ions& ions);

    void computeHvnlMatrix(
        const Ions&, dist_matrix::SparseDistMatrix<DISTMATDTYPE>&) const;
    void computeHvnlMatrix(const Ions&, ProjectedMatricesInterface*) const;

    void computeHvnlMatrix(const KBPsiMatrixInterface* const, const Ions&,
        dist_matrix::SparseDistMatrix<DISTMATDTYPE>&) const;
    void computeHvnlMatrix(const KBPsiMatrixInterface* const, const Ions&,
        ProjectedMatricesInterface*) const;

    void computeAll(Ions& ions, LocGridOrbitals& orbitals);
    void setup(const Ions& ions, const LocGridOrbitals& orbitals);
};

#endif
