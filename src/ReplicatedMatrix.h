#ifndef MGMOL_REPLICATEDMATRIX_H
#define MGMOL_REPLICATEDMATRIX_H

#ifdef HAVE_MAGMA

#include <string>
#include <vector>
#include <memory>

class ReplicatedMatrix
{
    // matrix size
    int dim_;

    // leading dimension for storage
    size_t ld_;

    // matrix data
    std::unique_ptr<double, void (*)(double*)> device_data_;

public:
    ReplicatedMatrix(const std::string name, const int m, const int n);
    ReplicatedMatrix(const std::string name, const int n);

    ReplicatedMatrix(const ReplicatedMatrix&);

    ~ReplicatedMatrix();

    ReplicatedMatrix& operator-=(const ReplicatedMatrix& rhs)
    {
        axpy(-1.0, rhs);
        return *this;
    }
    ReplicatedMatrix& operator=(const ReplicatedMatrix& rhs);

    void axpy(const double alpha, const ReplicatedMatrix& a);

    void setRandom(const double minv, const double maxv);
    void identity();
    void transpose(
        const double alpha, const ReplicatedMatrix&, const double beta);
    void trmm(const char, const char, const char, const char, const double,
        const ReplicatedMatrix&);
    void trtrs(const char, const char, const char, ReplicatedMatrix&) const;

    int potrf(char uplo);
    int potri(char uplo);
    void potrs(char, ReplicatedMatrix&);
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

    void print(
        std::ostream& os, const int, const int, const int, const int) const;
    void printMM(std::ostream& os) const;

    void clear(void);
    void trset(const char uplo);
};

void rotateSym(ReplicatedMatrix&, const ReplicatedMatrix&, ReplicatedMatrix&);

#endif

#endif
