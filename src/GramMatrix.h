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

#include "DistMatrix.h"
#include "MGmol_MPI.h"

#define NPRINT_ROWS_AND_COLS 5

class GramMatrix
{
    unsigned int dim_;

    dist_matrix::DistMatrix<DISTMATDTYPE>* matS_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* invS_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* ls_;

    dist_matrix::DistMatrix<DISTMATDTYPE>* work_;

    // index of orbitals associated to object (for compatibility testing)
    int orbitals_index_;
    bool isLSuptodate_;
    bool isInvSuptodate_;

    void transformLTML(dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
        const DISTMATDTYPE alpha) const;

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

    const dist_matrix::DistMatrix<DISTMATDTYPE>& getMatrix() const
    {
        assert(matS_ != nullptr);
        return *matS_;
    }

    const dist_matrix::DistMatrix<DISTMATDTYPE>& getInverse() const
    {
        assert(isInvSuptodate_);

        return *invS_;
    }

    const dist_matrix::DistMatrix<DISTMATDTYPE>& getCholeskyL() const
    {
        assert(ls_ != NULL);
        assert(isLSuptodate_);

        return *ls_;
    }

    void setMatrix(const dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
        const int orbitals_index);

    void updateLS();

    void set2Id(const int orbitals_index);

    void applyInv(dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
    {
        assert(isLSuptodate_);

        ls_->potrs('l', mat);
    }

    void printMM(std::ostream& tfile) { matS_->printMM(tfile); }

    void computeInverse();
    void solveLST(dist_matrix::DistMatrix<DISTMATDTYPE>& z) const;
    double computeCond();
    void sygst(dist_matrix::DistMatrix<DISTMATDTYPE>& mat);
    void rotateAll(const dist_matrix::DistMatrix<DISTMATDTYPE>& matU);

    double getLinDependent2states(int& st1, int& st2, int& st3) const;
    double getLinDependent2states(int& st1, int& st2) const;

    void computeLoewdinTransform(
        dist_matrix::DistMatrix<DISTMATDTYPE>& loewdinMat,
        std::shared_ptr<dist_matrix::DistMatrix<DISTMATDTYPE>> invLoewdin,
        const int orb_index);
};

#endif
