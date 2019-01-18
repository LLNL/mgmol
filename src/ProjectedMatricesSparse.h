// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_PROJECTED_MATRICES_SPARSE_H
#define MGMOL_PROJECTED_MATRICES_SPARSE_H

#include "ProjectedMatricesInterface.h"

#include <string.h>

#include "ClusterOrbitals.h"
#include "DensityMatrixSparse.h"
#include "LinearSolver.h"
#include "LocalizationRegions.h"
#include "MPIdata.h"
#include "ShortSightedInverse.h"
#include "SquareLocalMatrices.h"
#include "Timer.h"
#include "VariableSizeMatrix.h"

#include "hdf5.h"
#include <iostream>

// static const int MAX_MAT_SIZE = 10000;
// const double tol_matrix_elements=1.e-14;
// const double tol_matrix_elements=0.;

class HDFrestart;

// struct to hold a double and the associated rank
struct typeDoubleInt
{
    double val;
    int rank;
};

class ProjectedMatricesSparse : public ProjectedMatricesInterface
{
    bool isDataSetup_;

    double min_val_;

    int iterative_index_h_;
    int orbitals_index_;

    int index(const int phi_index, const int pot_index) const
    {
        return (phi_index % 100) + 100 * pot_index;
    }

    static Timer compute_inverse_tm_;
    static Timer compute_invB_tm_;
    static Timer init_gram_matrix_tm_;
    static Timer update_theta_tm_;
    static Timer update_submatT_tm_;
    static Timer update_submatX_tm_;
    static Timer eigsum_tm_;
    static Timer consolidate_H_tm_;
    static Timer eig_interval_tm_;

    ShortSightedInverse* invS_;
    VariableSizeMatrix<sparserow>* submatT_;
    VariableSizeMatrix<sparserow>* sH_;
    VariableSizeMatrix<sparserow>* matHB_;

    DensityMatrixSparse* dm_;

    SquareLocalMatrices<MATDTYPE>* localX_;
    SquareLocalMatrices<MATDTYPE>* localT_;

    /* Data distribution objects */
    DataDistribution* distributor_sH_;
    DataDistribution* distributor_matS_;
    DataDistribution* distributor_invS_;

    void clearData();

    int dim_;

    /* local data size */
    int lsize_;

    LocalizationRegions* lrs_;
    ClusterOrbitals* local_cluster_;

    // indexes corresponding to valid function in each subdomain
    std::vector<std::vector<int>> global_indexes_;

    std::vector<int> locvars_;

    void updateLocalMat(const VariableSizeMatrix<sparserow>& submatM,
        SquareLocalMatrices<MATDTYPE>* localM);
    double getExpectationMat(VariableSizeMatrix<sparserow>* mat);
    VariableSizeMatrix<sparserowtab>* getGramMat() { return invS_->gramMat(); }
    void computeGenEigenInterval(std::vector<double>& interval,
        const int maxits, const double padding = 0.01);

    double eigenvalue0_;

public:
    ProjectedMatricesSparse(const int ndim, LocalizationRegions* lrs,
        ClusterOrbitals* local_cluster = 0);
    ~ProjectedMatricesSparse();

    virtual void setup(const double kbt, const int nel,
        const std::vector<std::vector<int>>& global_indexes);

    void updateSubMatT();
    void updateTheta();
    double getExpectationH();
    void consolidateH();
    void consolidateOrbitalsOverlapMat(VariableSizeMatrix<sparserow>& mat);
    double dotProductWithInvS(const SquareLocalMatrices<MATDTYPE>& ss);
    double dotProductSimple(const SquareLocalMatrices<MATDTYPE>& ss);
    double dotProductWithDM(const SquareLocalMatrices<MATDTYPE>& ss);
    void applyInvS(SquareLocalMatrices<MATDTYPE>& mat);
    double getEigSum();
    double getLinDependent2states(
        int& st1, int& st2, const bool print_flag = false) const;
    void printMatrices(ostream& os) const;
    void printTimers(ostream& os);
    double getNel() const;
    double getLowestEigenvalue() const
    {
        assert(eigenvalue0_ == eigenvalue0_);
        return eigenvalue0_;
    }
    void printGramMatrix2states(
        const int st1, const int st2, std::ostream& os) const;

    void initializeGramMatrix(
        const SquareLocalMatrices<MATDTYPE>& ss, const int orbitals_index)
    {
        assert(invS_ != NULL);
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.barrier();
        init_gram_matrix_tm_.start();
        (*invS_).initGramMatrix(ss, global_indexes_, orbitals_index);

        init_gram_matrix_tm_.stop();
    }

    SquareLocalMatrices<MATDTYPE>& getLocalX() const { return *localX_; }

    SquareLocalMatrices<MATDTYPE>& getLocalT() const { return *localT_; }

    void updateSubMatX()
    {
        assert(dm_ != NULL);
        update_submatX_tm_.start();
        (*dm_).getLocalMatrix((*localX_), global_indexes_);
        update_submatX_tm_.stop();
    }

    void updateHB()
    {
        /* scale H */
        //    (*matHB_).scale(vel_);
    }
    void setDMuniform(const PROJMATDTYPE nel, const int orbitals_index)
    {
        dm_->setUniform(nel, orbitals_index);
    }

    int getDMMatrixIndex() const
    {
        assert(dm_ != 0);
        return (*dm_).getOrbitalsIndex();
    }
    int getGramMatrixIndex() const
    {
        assert(invS_ != 0);
        return (*invS_).getGramMatrixOrbitalsIndex();
    }
    double computeEntropy() { return 0.; }

    // This call assumes sH has not been consolidated and its row ordering on
    // the local subdomain (locvars_) has been defined by calling
    // setupSparseRows(locvars_). NOTE: We restrict only to contributions
    // computed/ belonging in the local subdomain
    void addMatrixElementSparseH(const int st1, const int st2, const double val)
    {
        assert(sH_->n() == (int)locvars_.size());

        (*sH_).insertMatrixElement(st1, st2, val, ADD, false);

        return;
    }

    // Add data from square local matrix (only contributions from functions
    // overlapping subdomain)
    void addMatrixElementsSparseH(const SquareLocalMatrices<MATDTYPE>& slH)
    {

        Control& ct = *(Control::instance());
        (*sH_).initializeMatrixElements(
            slH, global_indexes_, ct.numst, tol_matrix_elements);

        return;
    }

    void clearSparseH()
    {
        (*sH_).reset();
        /* initialize sparse H */
        (*sH_).setupSparseRows(locvars_);
        return;
    }

    void scaleSparseH(const double scale)
    {
        (*sH_).scale(scale);

        return;
    }

    double getDMEntry(const int row, const int col) const
    {
        assert(dm_ != 0);

        return (*dm_).getEntry(row, col);
    }

    void getDMEntries(const int row, const std::vector<int>& cols,
        std::vector<double>& values) const
    {
        assert(dm_ != 0);

        return dm_->getEntries(row, cols, values);
    }

    /* print the Theta matrix -- for diagnostics */
    void printTheta(ostream& os) const
    {
        if (onpe0) os << " Matrix Theta ... n = " << (*submatT_).n() << endl;

        std::vector<int> locfcns;
        (*lrs_).getLocalSubdomainIndices(locfcns);
        (*submatT_).print(os, locfcns);

        return;
    }

    /* print the HB matrix -- for diagnostics */
    void printHB(ostream& os) const
    {

        if (onpe0) os << " Matrix HB ... n = " << (*matHB_).n() << endl;

        std::vector<int> locfcns;
        (*lrs_).getLocalSubdomainIndices(locfcns);
        (*matHB_).print(os, locfcns);

        return;
    }

    /* print the inverse Gram matrix -- for diagnostics */
    void printDM(ostream& os) const
    {
        assert(dm_ != NULL);
        (*dm_).printDM(os);

        return;
    }

    void setDMto2InvS()
    {
        assert(invS_ != NULL);
        (*dm_).setto2InvS(
            (*invS_).getInvS(), (*invS_).getGramMatrixOrbitalsIndex());

        return;
    }

    void computeInvS()
    {
        assert(invS_ != NULL);
        assert(distributor_invS_ != NULL);
        assert(distributor_matS_ != NULL);

        compute_inverse_tm_.start();
        (*invS_).computeInvS(*distributor_matS_, *distributor_invS_);

        //   double rcond = checkCond(1000.);
        //   if(onpe0)cout<<"RCOND = "<<rcond<<endl;
        compute_inverse_tm_.stop();
    }
    void updateThetaAndHB()
    {
        updateTheta();
        updateHB();
        return;
    }

    void printGramMM(ofstream& tfile)
    {
        assert(invS_ != NULL);
        (*invS_).printGramMM(tfile);

        return;
    }

    /* print (part of) the (global) Gram matrix -- for diagnostics */
    void printS(ostream& os) const
    {
        assert(invS_ != NULL);
        (*invS_).printGramMat(os);
    }

    /* print the H matrix -- for diagnostics */
    void printH(ostream& os) const
    {
        /* matrix sH has local contributions, so we need to
         * gather neighboring data first.
         */
        if (onpe0) os << " Matrix H " << endl;

        std::vector<int> locfcns;
        (*lrs_).getLocalSubdomainIndices(locfcns);
        (*sH_).print(os, locfcns);

        return;
    }

    void resetDotProductMatrices()
    {
        assert((*invS_).gramMat() != 0);

        return;
    }

    void computeInvB()
    {
        assert(invS_ != NULL);
        compute_invB_tm_.start();
        (*invS_).computeInvS(*distributor_matS_, *distributor_invS_);
        compute_invB_tm_.stop();
        return;
    }

    double computeCond()
    {
        assert(invS_ != NULL);

        return (*invS_).computeGramMatCond();
    }

    double checkCond(const double tol, const bool flag = true)
    {
        double rcond = computeCond();

        if (rcond > tol)
        {
            // ofstream tfile("s.mm", ios::out);;
            // gm_->printMM(tfile);
            // tfile.close();
#ifdef USE_MPI
            MGmol_MPI& mgmolmpi = *(MGmol_MPI::instance());
            mgmolmpi.barrier();
#endif
            if (onpe0)
                (*MPIdata::sout)
                    << " CONDITION NUMBER OF THE OVERLAP MATRIX EXCEEDS TOL: "
                    << rcond << "!!!" << endl;
            Control& ct = *(Control::instance());
            if (flag) ct.global_exit(2);
        }
        return rcond;
    }

    double getTraceDiagProductWithInvS(std::vector<PROJMATDTYPE>& ddiag)
    {
        assert(invS_ != NULL);
        return (*invS_).getTraceDiagProductWithInvS(ddiag);
    }

    bool precondIsNew() const { return invS_->precondIsNew(); }

    bool isGramMatAugmented() { return invS_->isGramAugmented(); }

    double computeTraceInvSmultMat(
        const SquareLocalMatrices<MATDTYPE>& mat, const bool consolidate);
    double computeTraceDMmultMat(
        const SquareLocalMatrices<MATDTYPE>& mat, const bool consolidate);
    double computeTraceInvSmultMatmultTheta(
        SquareLocalMatrices<MATDTYPE>& mat, const bool consolidate);
    double computeTraceInvSmultMat(
        VariableSizeMatrix<sparserow>& vsmat, const bool consolidate);
    virtual void stripDM() { return; };
    virtual void dressupDM() { return; };

    const VariableSizeMatrix<sparserow>& getH() { return *matHB_; }
    void assignH(const VariableSizeMatrix<sparserow>& matH) { *matHB_ = matH; }
    void printEigenvalues(ostream& os) const {}
    void printOccupations(ostream& os) const {}
    void setHB2H() {}

    DensityMatrixSparse& getDM() { return *dm_; }
};

#endif
