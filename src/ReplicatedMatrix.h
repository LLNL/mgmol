#ifndef MGMOL_REPLICATEDMATRIX_H
#define MGMOL_REPLICATEDMATRIX_H

#ifdef HAVE_MAGMA

void rotateSym(ReplicatedMatrix&, const ReplicatedMatrix&, ReplicatedMatrix&);

class ReplicatedMatrix
{
    // matrix size
    int dim_;

    // leading dimension for storage
    size_t ld_;

    // matrix data
    double* device_data_;

public:
    ReplicatedMatrix(const std::string name, const int m, const int n);

    ReplicatedMatrix(const ReplicatedMatrix&);

    void identity();

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

    void print(
        std::ostream& os, const int, const int, const int, const int) const;
    void printMM(std::ostream& os) const;

    void clear(void);
};

#endif

#endif
