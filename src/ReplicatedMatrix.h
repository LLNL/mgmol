// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef MGMOL_REPLICATEDMATRIX_H
#define MGMOL_REPLICATEDMATRIX_H

#ifdef HAVE_MAGMA

class ReplicatedVector;
#include "SquareLocalMatrices.h"
#include "SquareSubMatrix.h"

#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

class ReplicatedMatrix
{
    static MPI_Comm comm_;
    static bool onpe0_;

    // matrix size
    int dim_;

    // leading dimension for storage
    size_t ld_;

    // matrix data
    std::unique_ptr<double, void (*)(double*)> device_data_;

    std::string name_;

public:
    friend class ReplicatedVector;

    static void setMPIcomm(MPI_Comm comm)
    {
        comm_ = comm;
        int mpi_rank;
        MPI_Comm_rank(comm_, &mpi_rank);
        if (mpi_rank == 0)
            onpe0_ = true;
        else
            onpe0_ = false;
    }

    ReplicatedMatrix(const std::string name, const int m, const int n);
    ReplicatedMatrix(const std::string name, const int n);

    // construct diagonal matrix from diagonal values
    ReplicatedMatrix(const std::string name, const double* const diagonal,
        const int m, const int n);

    ReplicatedMatrix(const ReplicatedMatrix&);

    ~ReplicatedMatrix();

    std::string name() { return name_; }

    double* const data() const { return device_data_.get(); }

    int m() const { return dim_; }

    int ld() const { return ld_; }

    void getsub(const ReplicatedMatrix& a, int m, int n, int ia, int ja);

    ReplicatedMatrix& operator-=(const ReplicatedMatrix& rhs)
    {
        axpy(-1.0, rhs);
        return *this;
    }
    ReplicatedMatrix& operator=(const ReplicatedMatrix& rhs);

    void assign(const ReplicatedMatrix& src, const int ib, const int jb);

    template <typename MemorySpaceType>
    void assign(SquareLocalMatrices<double, MemorySpaceType>& src);

    void add(const SquareSubMatrix<double>& src);

    // sum up values across MPI tasks
    void consolidate();

    void axpy(const double alpha, const ReplicatedMatrix& a);

    void init(const double* const ha, const int lda);
    void get(double* ha, const int lda) const;
    void getDiagonalValues(double* ha);
    void setRandom(const double minv, const double maxv);
    void identity();

    void scal(const double alpha);
    double trace(void) const;
    void transpose(
        const double alpha, const ReplicatedMatrix&, const double beta);
    void trmm(const char, const char, const char, const char, const double,
        const ReplicatedMatrix&);
    void trtrs(const char, const char, const char, ReplicatedMatrix&) const;

    int potrf(char uplo);
    int potri(char uplo);
    void potrs(char, ReplicatedMatrix&);
    void potrs(char, ReplicatedVector&);

    void getrf(std::vector<int>& ipiv);
    void getrs(char trans, ReplicatedMatrix& b, std::vector<int>& ipiv);

    void gemm(const char transa, const char transb, const double alpha,
        const ReplicatedMatrix& a, const ReplicatedMatrix& b,
        const double beta);
    void symm(const char side, const char uplo, const double alpha,
        const ReplicatedMatrix& a, const ReplicatedMatrix& b,
        const double beta);
    void syev(char, char, std::vector<double>&, ReplicatedMatrix&);
    void sygst(int, char, const ReplicatedMatrix&);

    void setVal(const int i, const int j, const double val);
    void setDiagonal(const std::vector<double>& diag_values);
    int iamax(const int j, double& val);
    double norm(char ty);
    double traceProduct(const ReplicatedMatrix&) const;

    void print(
        std::ostream& os, const int, const int, const int, const int) const;
    void printMM(std::ostream& os) const;

    void clear(void);
    void trset(const char uplo);
};

void rotateSym(ReplicatedMatrix&, const ReplicatedMatrix&, ReplicatedMatrix&);

#endif

#endif
