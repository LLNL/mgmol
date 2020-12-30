// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_SP2_H
#define MGMOL_SP2_H

#include "stdlib.h"
#include <iostream>

#ifdef HAVE_BML
extern "C"
{
#include <bml.h>
}
#endif

#include "SquareLocalMatrices.h"
#include "Timer.h"

class SP2
{

    // we actually store X*S in these matrices
#ifdef HAVE_BML
    bml_matrix_t* Xi_;
    bml_matrix_t* Xi_sq_;
#else
    SquareLocalMatrices<MATDTYPE, MemorySpace::Host>* Xi_;
    SquareLocalMatrices<MATDTYPE, MemorySpace::Host>* Xi_sq_;
#endif

    const double tol_;

    // distributed computation of trace (vs. serial/replicated)
    const bool distributed_;

    //  double trace1_; //Trace[Xi_]
    //  double trace2_; //Trace[Xi_sq_]
    double trace_[2];

    // Calculate A for the current Xi_
    int calcA(const int nel);
    // Iterate solver starting from Xi_
    void iterate(int A);

    std::vector<int> loc_ids_;

    static Timer getdm_tm_;

    void reduceSumTrace();

public:
    // Create SP2 object with a variable size matrix theta
    // ratio=R/(R_s)<=1, tol the tolerance of the solver
    SP2(const double tol, const bool distributed);
    ~SP2();

    // Calculate density matrix using SP2
    void solve(const int nel, const bool verbose);

    // Map variable size matrix to square matrix
    template <class T>
    void initializeLocalMat(const T& submatM, const double emin,
        const double emax, const std::vector<int>& loc_ids);

    // Map local dm (internal) into global one, and multiply by inverse(S)
    template <class T>
    void getDM(T& submatM, const T& invS);

    static void printTimers(std::ostream& os) { getdm_tm_.print(os); }
};
#endif
