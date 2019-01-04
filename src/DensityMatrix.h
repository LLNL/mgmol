// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DENSITYMATRIX_H
#define MGMOL_DENSITYMATRIX_H

#include "DistMatrix.h"
#include "MPIdata.h"

#include <ostream>
#include <vector>

#define DM_NPRINT_ROWS_AND_COLS 5

class DensityMatrix
{
    int dim_;
    std::vector<DISTMATDTYPE> occupation_;

    dist_matrix::DistMatrix<DISTMATDTYPE>* dm_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* kernel4dot_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* work_;

    int orbitals_index_;

    bool occ_uptodate_;
    bool uniform_occ_;
    bool stripped_;
    DISTMATDTYPE orbital_occupation_;

    DensityMatrix();

public:
    DensityMatrix(const int ndim);
    DensityMatrix(const DensityMatrix&);

    ~DensityMatrix();

    DensityMatrix& operator=(const DensityMatrix&);
    void setUniform(const DISTMATDTYPE nel, const int new_orbitals_index);

    int getOrbitalsIndex() const { return orbitals_index_; }

    bool occupationsUptodate() const { return occ_uptodate_; }
    bool fromUniformOccupations() const { return uniform_occ_; }
    DISTMATDTYPE getVal(const int i) const
    {
        assert(!stripped_);
        return dm_->val(i);
    }

    double dot(const dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
    {
        assert(!stripped_);
        return dm_->dot(mat);
    }

    void print(std::ostream& os) const
    {
        assert(!stripped_);
        if (onpe0) os << " DensityMatrix" << std::endl;
        dm_->print(os, 0, 0, DM_NPRINT_ROWS_AND_COLS, DM_NPRINT_ROWS_AND_COLS);
    }

    const dist_matrix::DistMatrix<DISTMATDTYPE>& getMatrix() const
    {
        assert(!stripped_);
        assert(dm_ != 0);
        return *dm_;
    }

    const dist_matrix::DistMatrix<DISTMATDTYPE>& kernel4dot() const
    {
        return *kernel4dot_;
    }

    void setMatrix(const dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
        const int orbitals_index)
    {
        *dm_            = mat;
        orbitals_index_ = orbitals_index;

        occupation_.clear();

        occ_uptodate_ = false;
        uniform_occ_  = false;
        stripped_     = false;
    }
    void initMatrix(const DISTMATDTYPE* const val)
    {
        dm_->init(val, dim_);
        occupation_.clear();

        occ_uptodate_ = false;
        uniform_occ_  = false;
        stripped_     = false;

        orbitals_index_ = 0;
    }
    // double getSumOccupations()const;
    void getOccupations(std::vector<DISTMATDTYPE>& occ) const
    {
        assert(occ_uptodate_);
        assert((int)occ.size() == dim_);
        memcpy(&occ[0], &occupation_[0], dim_ * sizeof(DISTMATDTYPE));
    }

    void setOccupations(const std::vector<DISTMATDTYPE>& occ);

    void setto2InvS(const dist_matrix::DistMatrix<DISTMATDTYPE>& invS,
        const int orbitals_index);

    void stripS(const dist_matrix::DistMatrix<DISTMATDTYPE>& ls,
        const int orbitals_index_gram);
    void dressUpS(const dist_matrix::DistMatrix<DISTMATDTYPE>& ls,
        const int new_orbitals_index);

    void buildFromBlock(const dist_matrix::DistMatrix<DISTMATDTYPE>& block00);

    double computeEntropy() const;
    void computeOccupations(const dist_matrix::DistMatrix<DISTMATDTYPE>& ls);
    void build(const int new_orbitals_index);
    void build(
        const std::vector<DISTMATDTYPE>& occ, const int new_orbitals_index);
    void build(const dist_matrix::DistMatrix<DISTMATDTYPE>& z,
        const int new_orbitals_index);
    void build(const dist_matrix::DistMatrix<DISTMATDTYPE>& z,
        const std::vector<DISTMATDTYPE>& occ, const int new_orbitals_index);

    void rotate(const dist_matrix::DistMatrix<DISTMATDTYPE>& rotation_matrix,
        const bool flag_eigen);
    void printOccupations(std::ostream& os) const;
    void diagonalize(const dist_matrix::DistMatrix<DISTMATDTYPE>& ls,
        std::vector<DISTMATDTYPE>& occ);
    void diagonalize(const char eigv, std::vector<DISTMATDTYPE>& occ,
        dist_matrix::DistMatrix<DISTMATDTYPE>& vect);
    double getExpectation(const dist_matrix::DistMatrix<DISTMATDTYPE>& A);
};

#endif
