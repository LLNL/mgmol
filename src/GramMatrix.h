// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_GRAMMATRIX_H
#define MGMOL_GRAMMATRIX_H

#include "MGmol_MPI.h"

#include <memory>

#define NPRINT_ROWS_AND_COLS 5

template <class MatrixType>
class GramMatrix
{
    unsigned int dim_;

    MatrixType* matS_;
    MatrixType* invS_;
    MatrixType* ls_;

    MatrixType* work_;

    // index of orbitals associated to object (for compatibility testing)
    int orbitals_index_;
    bool isLSuptodate_;
    bool isInvSuptodate_;

    void transformLTML(MatrixType& mat, const double alpha) const;

public:
    GramMatrix(const int ndim);
    GramMatrix(const GramMatrix&);
    GramMatrix& operator=(const GramMatrix&);

    ~GramMatrix();

    int getAssociatedOrbitalsIndex() const { return orbitals_index_; }

    void print(std::ostream& os, int nprint_rows = NPRINT_ROWS_AND_COLS,
        int nprint_col = NPRINT_ROWS_AND_COLS) const
    {
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        if (mmpi.instancePE0()) os << " GramMatrix" << std::endl;
        matS_->print(os, 0, 0, nprint_rows, nprint_col);
    }

    void printMM(std::ostream& os) const
    {
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        if (mmpi.instancePE0()) os << "Gram Matrix" << std::endl;
        matS_->printMM(os);
    }

    const MatrixType& getMatrix() const
    {
        assert(matS_ != nullptr);
        return *matS_;
    }

    const MatrixType& getInverse() const
    {
        assert(isInvSuptodate_);

        return *invS_;
    }

    const MatrixType& getCholeskyL() const
    {
        assert(ls_ != NULL);
        assert(isLSuptodate_);

        return *ls_;
    }

    void setMatrix(const MatrixType& mat, const int orbitals_index);

    void updateLS();

    void set2Id(const int orbitals_index);

    void applyInv(MatrixType& mat);
    template <class VectorType>
    void applyInv(VectorType& v);

    void printMM(std::ostream& tfile) { matS_->printMM(tfile); }

    void computeInverse();
    void solveLST(MatrixType& z) const;
    double computeCond();
    void sygst(MatrixType& mat);
    void rotateAll(const MatrixType& matU);

    double getLinDependent2states(int& st1, int& st2, int& st3) const;
    double getLinDependent2states(int& st1, int& st2) const;

    void computeLoewdinTransform(MatrixType& loewdinMat,
        std::shared_ptr<MatrixType> invLoewdin, const int orb_index);

    double getTraceDiagProductWithInvS(std::vector<double>& ddiag);
};

#endif
