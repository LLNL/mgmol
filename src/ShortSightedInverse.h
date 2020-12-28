// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
 * Main header file for inverting the Gram Matrix
 */

#ifndef _SHORTSIGHTEDINVERSE_H_
#define _SHORTSIGHTEDINVERSE_H_

#include <vector>

#include <mpi.h>

#include "ClusterOrbitals.h"
#include "DataDistribution.h"
#include "LinearSolver.h"
#include "LinearSolverMatrix.h"
#include "LocalizationRegions.h"
#include "Mesh.h"
#include "PreconILU.h"
#include "SquareLocalMatrices.h"
#include "VariableSizeMatrix.h"

typedef double DMDTYPE;

class ShortSightedInverse
{
    static Timer Gram_Matrix_data_distribution_tm_;
    static Timer reset_tm_;
    static Timer inverse_solve_tm_;
    static Timer compute_invS_tm_;
    static Timer gather_invS_tm_;
    static Timer gmres_tm_;
    static Timer pc_setup_tm_;
    static Timer pc_solve_tm_;
    static Timer linear_solver_matvec_tm_;
    static Timer linear_solver_matrix_init_tm_;

    /* Local variables */
    int gm_orbitals_index_;
    bool issetup_; // flag for initial data setup
    bool isInvSUpToDate_; // flag for checking if inverse is up to date
    bool isGramAugmented_; // flag for checking if gram matrix is augmented
    const std::vector<int>&
        locvars_; // array of global indexes overlapping with this PE

    /* Matrix variables */
    int lsize_; // Initial size of local matrix
    int aug_size_; // Augmented size of the local matrix after data distribution
    int max_size_; // max. local matrix size
    bool recompute_pc_; // flag to recompute the preconditioner (or not)
    bool new_pc_; // preconditioner recompute in last solve
    VariableSizeMatrix<sparserowtab>* gramMat_; // Matrix for data distribution
    VariableSizeMatrix<sparserow>*
        invS_; // Matrix for storing gram matrix inverse
    LinearSolverMatrix<lsdatatype>*
        matLS_; // Linear solver matrix for linear solver
    PreconILU<pcdatatype>* precon_; // preconditioner for linear system solve
    /* Preconditioner parameters */
    float droptol_;
    short MaxFil_;
    short lof_;
    PCILUTYPE ilutype_;

    /* Data Distribution variables */
    double loc_radius_; // Localization radius
    //  double spread_radius_;				// Spreading radius for
    //  data distribution

    /* Linear system solver variables */
    std::vector<int>
        locfcns_; // Global ids of Local variables centered on this proc.
    std::vector<int>
        loccluster_; // indexes of cluster assigned to this subdomain
    double fgmres_tol_; // relative convergence tolerance for fgmres
    short im_; // Krylov subspace dimension for fgmres
    short maxits_; // max. fgmres iterations allowed
    double resnorm_; // fgmres residual norm

    void gather(DataDistribution& dtor); // gather inverse data
    void reset(); // Reset/ setup data

    void augmentGramMatrix();

    //  int solve(LinearSolverMatrix<lsdatatype>& LSMat);
    int solve();

public:
    ShortSightedInverse(std::shared_ptr<LocalizationRegions> lrs,
        const std::vector<int>& locvars,
        ClusterOrbitals* local_cluster = nullptr); // constructor
    void initGramMatrix(const LocalMatrices<MATDTYPE, MemorySpace::Host>& ss,
        const std::vector<std::vector<int>>& global_indexes,
        const int orbitals_index); // initialize local Gram matrix
    void computeInvS(
        DataDistribution& dtor_invS); // Compute the local matrix inverse
    double getTraceDiagProductWithInvS(
        std::vector<DMDTYPE>& ddiag); // return sum_i ( ddiag[i]*S**(-1)[i][i] )
    void printGramMat(std::ostream& os, int nrows = NUM_PRINT_ROWS) const;
    double getTraceDotProductWithInvS(VariableSizeMatrix<sparserow>* mat);
    static void printTimers(std::ostream& os); // print timers
    ~ShortSightedInverse(); // destructor

    /* get the initial local size */
    int get_lsize() { return lsize_; }

    int getGramMatrixOrbitalsIndex() { return gm_orbitals_index_; }
    // matrix multiplication S**(-1) * B
    // flag== true => compute entries for specific nonzero pattern only
    void invSmultB(VariableSizeMatrix<SparseRow>* B,
        VariableSizeMatrix<SparseRow>& C, bool flag = true)
    {
        assert(B != NULL);

        (*invS_).AmultSymBLocal(B, C, locfcns_, *gramMat_, flag);

        return;
    }

    // return [S**(-1) * B]_ij
    double invSmultB_ij(
        VariableSizeMatrix<SparseRow>* B, const int row, const int col)
    {
        assert(invS_ != 0);
        assert((*invS_).n() != 0);

        return (*invS_).AmultSymB_ij(B, row, col);
    }

    void printGramMM(std::ofstream& tfile);

    /* compute condition number of gramMat_ */
    double computeGramMatCond();

    VariableSizeMatrix<sparserow>& getInvS() { return *invS_; }

    bool precondIsNew() const { return new_pc_; }

    bool isGramAugmented() const { return isGramAugmented_; }

    bool isInvSUpToDate() { return isInvSUpToDate_; }

    VariableSizeMatrix<sparserowtab>* gramMat() { return gramMat_; }

    VariableSizeMatrix<sparserowtab>& getGramMat() { return *gramMat_; }

    std::vector<int> centeredFcnLocalIds();

    // Do local solve with (augmented) Gram Matrix
    int GramMatLSSolve(const double* rhs, double* sol);

    int augMatSize() { return aug_size_; }

    // get pointer to row of invS.
    double* getInvSRowEntries(const int row)
    {
        assert(row <= invS_->n());

        return invS_->getRowEntries(row);
    }

    PreconILU<pcdatatype>* precon() { return precon_; }

    LinearSolverMatrix<lsdatatype>* matLS() { return matLS_; }
};

#endif
