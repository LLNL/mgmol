// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_KBPSIMATRIX_SPARSE_H
#define MGMOL_KBPSIMATRIX_SPARSE_H

#include "DataDistribution.h"
#include "DensityMatrixSparse.h"
#include "DistMatrixWithSparseComponent.h"
#include "KBPsiMatrixInterface.h"
#include "VariableSizeMatrix.h"

#include "GridFunc.h"
#include "Lap.h"
#include "Timer.h"

#include <cassert>
#include <limits.h>

class Ions;
class Ion;
class LocGridOrbitals;
class ProjectedMatricesInterface;
class ProjectedMatricesSparse;

class KBPsiMatrixSparse : public KBPsiMatrixInterface
{
    static Timer global_sum_tm_;
    static Timer compute_kbpsi_tm_;
    static Timer computeHvnlMatrix_tm_;
    static Timer setup_tm_;
    static Timer trace_tm_;

    // operator used to compute r.h.s
    // usually identity, unless Mehrstellen scheme is used
    const pb::Lap<ORBDTYPE>* lapop_;

    const bool need2radius_;

    bool isDataSetup_;

    // Sparse matrix - VariableSizeMatrix
    VariableSizeMatrix<sparserow>* kbpsimat_; // matrix <KB, psi>
    VariableSizeMatrix<sparserow>* kbBpsimat_; // matrix <KB, B*psi>
    DataDistribution* distributor_; // data distribution object
    double spread_radius_; // radius for spreading data

    int numst_; // global number of functions

    int count_proj_subdomain_;

    void addKBPsi(const int gid, const int st, const double val)
    {
        assert(st < numst_);
        (*kbpsimat_).insertMatrixElement(gid, st, val, ADD, true);
    }
    void addKBBPsi(const int gid, const int st, const double val)
    {
        assert(st < numst_);
        (*kbBpsimat_).insertMatrixElement(gid, st, val, ADD, true);
    }
    double getKBPsi(const int gid, const int st) const
    {
        return (*kbpsimat_).get_value(gid, st);
    }
    double getKBBPsi(const int gid, const int st) const
    {
        return (*kbBpsimat_).get_value(gid, st);
    }

    void reset()
    {
        if (kbpsimat_ != 0) (*kbpsimat_).reset();

        if (lapop_)
            if (kbBpsimat_ != 0) (*kbBpsimat_).reset();
    }

    void globalSumKBpsi();

    void computeHvnlMatrix(const KBPsiMatrixSparse* const kbpsi, const Ion&,
        dist_matrix::SparseDistMatrix<DISTMATDTYPE>&) const;
    void computeHvnlMatrix(const KBPsiMatrixSparse* const kbpsi2,
        const Ion& ion, VariableSizeMatrix<sparserow>& mat) const;
    void computeHvnlMatrix(const KBPsiMatrixSparse* const kbpsi2, const Ion&,
        ProjectedMatricesInterface*) const;

    void getPsiKBPsiSym(const Ions& ions, VariableSizeMatrix<sparserow>& sm);
    void computeKBpsi(Ions& ions, LocGridOrbitals& orbitals,
        const int first_color, const int nb_colors, const bool flag);
    void clearData();

public:
    KBPsiMatrixSparse(pb::Lap<ORBDTYPE>* lapop, const bool need2radius = true);

    ~KBPsiMatrixSparse();

    void printTimers(ostream& os);

    double getEvnl(const Ions& ions, LocGridOrbitals& orbitals,
        ProjectedMatricesSparse* proj_matrices);
    void computeKBpsi(
        Ions& ions, pb::GridFunc<ORBDTYPE>*, const int, const bool flag);
    double getValIonState(const int gid, const int st) const
    {
        assert(st < numst_);
        return (*kbpsimat_).get_value(gid, st);
    }
    void scaleWithKBcoeff(const Ions& ions);

    void computeHvnlMatrix(
        const Ions&, dist_matrix::SparseDistMatrix<DISTMATDTYPE>&) const;
    void computeHvnlMatrix(const Ions&, VariableSizeMatrix<sparserow>&) const;
    void computeHvnlMatrix(const Ions&, ProjectedMatricesInterface*) const;

    void computeHvnlMatrix(const KBPsiMatrixInterface* const kbpsi, const Ions&,
        dist_matrix::SparseDistMatrix<DISTMATDTYPE>&) const;
    void computeHvnlMatrix(const KBPsiMatrixInterface* const kbpsi2,
        const Ions& ions,
        dist_matrix::DistMatrixWithSparseComponent<DISTMATDTYPE>& Aij) const;
    void computeHvnlMatrix(const KBPsiMatrixInterface* const kbpsi, const Ions&,
        VariableSizeMatrix<sparserow>&) const;
    void computeHvnlMatrix(const KBPsiMatrixInterface* const kbpsi2,
        const Ions& ions, ProjectedMatricesInterface* proj_matrices) const;

    void computeAll(Ions& ions, LocGridOrbitals& orbitals);
    void setup(const Ions& ions, const LocGridOrbitals& orbitals);

    double getTraceDM(
        const int gid, const DISTMATDTYPE* const mat_X, const int numst) const;
    double getTraceDM(const int gid, const DensityMatrixSparse& dm) const;
};

#endif
