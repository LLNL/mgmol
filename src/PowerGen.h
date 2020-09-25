// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_POWERGEN_H
#define MGMOL_POWERGEN_H

#include "GramMatrix.h"
#include "random.h"

#include <vector>

template <class MatrixType, class VectorType>
class PowerGen
{
private:
    static Timer compute_tm_;

    // use shift to target highest or lowest eigenvalue
    double shift_;

    VectorType vec1_;
    VectorType vec2_;

public:
    PowerGen(const int n)
        : shift_(0.), vec1_(generate_rand(n)), vec2_(generate_rand(n))
    {
    }

    void computeGenEigenInterval(MatrixType& mat, GramMatrix<MatrixType>& gm,
        std::vector<double>& interval, const int maxits, const double pad);

    static void printTimers(std::ostream& os) { compute_tm_.print(os); }
};

template <class MatrixType, class VectorType>
Timer PowerGen<MatrixType, VectorType>::compute_tm_("PowerGen::compute");
#endif
