// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_KBPSIMATRIX_SPARSE_H
#define MGMOL_KBPSIMATRIX_SPARSE_H

#include "DataDistribution.h"
#include "DensityMatrixSparse.h"
#include "Ions.h"
#include "KBPsiMatrixInterface.h"
#include "ProjectedMatricesSparse.h"
#include "SquareSubMatrix.h"
#include "VariableSizeMatrix.h"

#include "GridFunc.h"
#include "Lap.h"
#include "Timer.h"

#include <cassert>
#include <limits.h>

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

    int count_proj_subdomain_;

    void addKBPsi(const int gid, const int st, const double val) override
    {
        (*kbpsimat_).insertMatrixElement(gid, st, val, ADD, true);
    }
    void addKBBPsi(const int gid, const int st, const double val) override
    {
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

    // private functions working with single Ion
    void computeHvnlElementsIon(const KBPsiMatrixSparse* const kbpsi,
        const Ion&, SquareSubMatrix<double>& mat) const;
    void computeHvnlElementsIon(const KBPsiMatrixSparse* const kbpsi2,
        const Ion& ion, VariableSizeMatrix<sparserow>& mat) const;

    template <class OrbitalsType>
    void computeKBpsi(const Ions& ions, OrbitalsType& orbitals,
        const int first_color, const int nb_colors, const bool flag);
    void clearData();

public:
    KBPsiMatrixSparse(pb::Lap<ORBDTYPE>* lapop, const bool need2radius = true);

    ~KBPsiMatrixSparse() override;

    void printTimers(std::ostream& os) override;

    void reset()
    {
        if (kbpsimat_ != nullptr) (*kbpsimat_).reset();

        if (lapop_)
            if (kbBpsimat_ != nullptr) (*kbBpsimat_).reset();
    }

    void globalSumKBpsi();

    double getEvnl(const Ions& ions, ProjectedMatricesSparse* proj_matrices);
    template <class MatrixType>
    double getEvnl(
        const Ions& ions, ProjectedMatrices<MatrixType>* proj_matrices);
    void computeKBpsi(
        const Ions& ions, pb::GridFunc<ORBDTYPE>*, const int, const bool flag);
    double getValIonState(const int gid, const int st) const
    {
        return (*kbpsimat_).get_value(gid, st);
    }
    void scaleWithKBcoeff(const Ions& ions);

    SquareSubMatrix<double> computeHvnlMatrix(const Ions&) const;
    void computeHvnlMatrix(const Ions&, VariableSizeMatrix<sparserow>&) const;
    void computeHvnlMatrix(const Ions&, ProjectedMatricesInterface*) const;

    SquareSubMatrix<double> computeHvnlMatrix(
        const KBPsiMatrixInterface* const kbpsi, const Ions&) const;
    template <class MatrixType>
    void computeHvnlMatrix(const KBPsiMatrixInterface* const kbpsi2,
        const Ions& ions, MatrixType& Aij) const;
    void computeHvnlMatrix(const KBPsiMatrixInterface* const kbpsi, const Ions&,
        VariableSizeMatrix<sparserow>&) const;
    void computeHvnlMatrix(const KBPsiMatrixInterface* const kbpsi2,
        const Ions& ions, ProjectedMatricesInterface* proj_matrices) const;

    template <class T>
    void computeAll(const Ions& ions, T& orbitals);
    void setup(const Ions& ions);

    double getTraceDM(
        const int gid, const double* const mat_X, const int numst) const;
    double getTraceDM(const int gid, const DensityMatrixSparse& dm) const;
};

#endif
