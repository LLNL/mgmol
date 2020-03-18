// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
 * ILU preconditioning for Krylov Solver
 */

#ifndef _PRECONILU_H_
#define _PRECONILU_H_

#include "LinearSolverMatrix.h"
#include "mputils.h"
#include <vector>
/* define different preconditioner options */
typedef enum PCILUTYPE
{
    PCILUK,
    PCMILUT,
    PCILUT,
    PCILU0
} PCILUTYPE;
/* define matrix data type */
typedef float pcdatatype;
/* define default lof */
const int level_of_fill = 0;

template <class T>
class PreconILU
{
    typedef typename std::vector<T>::iterator TvecIterator;
    typedef typename std::vector<T>::const_iterator const_TvecIterator;
    typedef std::vector<lsdatatype> LSvec;
    typedef LSvec::iterator LSvecIterator;

    static Timer pcilu_setup_tm_;
    static Timer pcilu_solve_tm_;

    int n_; // matrix size
    LinearSolverMatrix<T>* L_; // L part elements
    std::vector<T> D_; // diagonal elements
    LinearSolverMatrix<T>* U_; // U part elements
    T droptol_; // defines drop tolerance for ILU
    int MaxFil_; // Max fill allowed for each row/col of L and U
    int lof_; // Level of fill for iluk only
    PCILUTYPE pctype_; // type of preconditioner

    /* set parameters for ILU factorization */
    void setParams(const float droptol, const int maxfill, const int lof);
    /* split an array into two sorted parts */
    void qsplitC(T* a, int* ind, const int n, const int Ncut);
    void qsplitC(
        std::vector<T>& a, std::vector<int>& ind, const int n, const int Ncut);
    /* Symbolic factorization for ILU */
    int lofC(LinearSolverMatrix<lsdatatype>& csmat);
    /* ILUK factorization */
    int ilukC(LinearSolverMatrix<lsdatatype>& csmat);
    /* ILU0 factorization */
    template <typename T2>
    int ilu0(LinearSolverMatrix<T2>& csmat);
    //  template<typename T2>
    //  int ilu02( LinearSolverMatrix<T2>& csmat);
    /* Modified ILUT factorization */
    int milut(LinearSolverMatrix<lsdatatype>& csmat);
    /* Standard ILUT factorization */
    int ilut(LinearSolverMatrix<lsdatatype>& csmat);
    /* Diagonal (row norm) Preconditioning */
    int diag(LinearSolverMatrix<lsdatatype>& csmat);

public:
    PreconILU(LinearSolverMatrix<lsdatatype>& csmat, const float droptol,
        const int maxfill,
        const int lof = level_of_fill); // construct ILU struct
    void setup(LinearSolverMatrix<lsdatatype>& csmat,
        const PCILUTYPE type); // perform ILU factorization
    int nnz_ilu(); // number of nonzero entries in ILU precon
    void LUsolve(double* const y,
        double* const x) const; // triangular solve with L and U factors
    void LUsolve(float* const y,
        float* const x) const; // triangular solve with L and U factors
    static void printTimers(std::ostream& os); // print timers
    ~PreconILU(); // destructor

    void diagSolve(double* const y_, double* const x_)
    {
        /*----------------------------------------------------------------------
         *    performs a diagonal solve (diagonal matrix product).
         *    y  = right-hand-side
         *    x  = solution on return
         *--------------------------------------------------------------------*/

        /* -------------begin------------*/

        const int n = n_;
        for (int i = 0; i < n; i++)
            x_[i] = D_[i] * y_[i];
    }
    int n() { return n_; }
};

#endif
