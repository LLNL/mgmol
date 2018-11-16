// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_PROJECTED_MATRICES_H
#define MGMOL_PROJECTED_MATRICES_H

#include "DensityMatrix.h"
#include "GramMatrix.h"
#include "MPIdata.h"
#include "ProjectedMatricesInterface.h"
#include "SquareLocalMatrices.h"
#include "Timer.h"
#include "tools.h"

#include "DistMatrix.h"
#include "RemoteTasksDistMatrix.h"
#include "SparseDistMatrix.h"

#include "hdf5.h"

#include <iostream>

#define EXTRAPOLATE_H 1

// const double tol_matrix_elements=1.e-14;
// const double tol_matrix_elements=0.;

class HDFrestart;

class ProjectedMatrices : public ProjectedMatricesInterface
{
    static short n_instances_;

    static short sparse_distmatrix_nb_tasks_per_partitions_;
    static GramMatrix* gram_4dotProducts_;
    static DensityMatrix* dm_4dot_product_;

    // spin: 0 for ignoring spin, 1 for calculation with spin
    const bool with_spin_;

    std::vector<DISTMATDTYPE> eigenvalues_;

    /*!
     * matrices to save old values and enable mixing or reset
     */
    dist_matrix::DistMatrix<DISTMATDTYPE>* mat_X_old_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* mat_L_old_;

    // orthogonal tranformation
    dist_matrix::DistMatrix<DISTMATDTYPE>* u_;
#ifdef PROCRUSTES
    dist_matrix::DistMatrix<DISTMATDTYPE>* u0_;
    // matrix of permutation to apply for continuity with u0_
    dist_matrix::DistMatrix<DISTMATDTYPE>* p_;
#endif
    double min_val_;

    static Timer sygv_tm_;
    static Timer compute_inverse_tm_;
    static Timer compute_invB_tm_;
    static Timer init_gram_matrix_tm_;
    static Timer update_theta_tm_;
    static Timer update_submatT_tm_;
    static Timer update_submatX_tm_;
    static Timer eigsum_tm_;
    static Timer consolidate_H_tm_;

    static dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>*
        remote_tasks_DistMatrix_;

#if EXTRAPOLATE_H
    dist_matrix::DistMatrix<DISTMATDTYPE>* h_minus1_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* h_minus2_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* new_h_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* tmp_h_minus1_;
#endif

    ProjectedMatrices& operator=(const ProjectedMatrices& src);
    ProjectedMatrices(const ProjectedMatrices& pm);

#ifdef USE_DIS_MAT
    dist_matrix::SubMatricesIndexing<DISTMATDTYPE>* submat_indexing_;

    dist_matrix::SubMatrices<DISTMATDTYPE>* submatLS_;
    dist_matrix::SubMatrices<DISTMATDTYPE>* submatWork_;
#endif

    /*!
     * Matrices to be used to multiply orbitals on the right
     * (dimension equal to number of local colors)
     */
    SquareLocalMatrices<MATDTYPE>* localX_; // density matrix
    SquareLocalMatrices<MATDTYPE>* localT_; // theta=inv(S)*H_phi

    dist_matrix::SparseDistMatrix<DISTMATDTYPE>* sH_;

    void printEigenvaluesHa(std::ostream& os) const;
    void printEigenvaluesEV(std::ostream& os) const;

    double computeChemicalPotentialAndOccupations(
        const std::vector<DISTMATDTYPE>& energies, const double width,
        const int nel, const int max_numst);

    double computeChemicalPotentialAndOccupations()
    {
        return computeChemicalPotentialAndOccupations(width_, nel_, dim_);
    }

protected:
    // indexes corresponding to valid function in each subdomain
    std::vector<std::vector<int>> global_indexes_;

    unsigned dim_;

    dist_matrix::DistMatrix<DISTMATDTYPE>* theta_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* matHB_;
    dist_matrix::DistMatrix<DISTMATDTYPE>* matH_;

    DensityMatrix* dm_;
    GramMatrix* gm_;

    // work matrix for tmp usage
    dist_matrix::DistMatrix<DISTMATDTYPE>* work_;

    std::vector<DISTMATDTYPE> aux_energies_;
    void printTheta(ostream& os) const
    {
        if (onpe0) os << " Matrix Theta" << endl;
        theta_->print(os, 0, 0, NPRINT_ROWS_AND_COLS, NPRINT_ROWS_AND_COLS);
    }

    void printHB(ostream& os) const
    {
        if (onpe0) os << " Matrix HB" << endl;
        matHB_->print(os, 0, 0, NPRINT_ROWS_AND_COLS, NPRINT_ROWS_AND_COLS);
    }

    void printH(ostream& os) const
    {
        if (onpe0) os << " Matrix H" << endl;
        matH_->print(os, 0, 0, NPRINT_ROWS_AND_COLS, NPRINT_ROWS_AND_COLS);
    }
    void setAuxilliaryEnergiesFromEigenenergies();

public:
    ProjectedMatrices(const int, const bool with_spin);
    virtual ~ProjectedMatrices();

    void setup(const double kbt, const int nel,
        const std::vector<std::vector<int>>& global_indexes);

    void addMatrixElementSparseH(const int st1, const int st2, const double val)
    {
        sH_->push_back(st1, st2, val);
    }

    // fill SparseDistMatrix sH_ with values in slH
    void addMatrixElementsSparseH(const SquareLocalMatrices<MATDTYPE>& slH)
    {
        Control& ct = *(Control::instance());
        slH.fillSparseDistMatrix(*sH_, global_indexes_, ct.numst);
    }
    void clearSparseH() { sH_->clearData(); }
    void scaleSparseH(const double scale) { sH_->scal(scale); }

    void consolidateH()
    {
        consolidate_H_tm_.start();
        sH_->parallelSumToDistMatrix();
        consolidate_H_tm_.stop();
    }

    double getTraceDiagProductWithInvS(std::vector<DISTMATDTYPE>& ddiag);

    SquareLocalMatrices<MATDTYPE>& getLocalX() const { return *localX_; }

    SquareLocalMatrices<MATDTYPE>& getLocalT() const { return *localT_; }

    const dist_matrix::SubMatrices<DISTMATDTYPE>& getSubMatLS() const
    {
        assert(submatLS_ != 0);
        return *submatLS_;
    }

    void updateSubMatX() { updateSubMatX(dm_->getMatrix()); }

    void updateSubMatX(const dist_matrix::DistMatrix<DISTMATDTYPE>& dm)
    {
        update_submatX_tm_.start();
        submatWork_->gather(dm);

        if (chromatic_number_ > 0)
            for (short iloc = 0; iloc < subdiv_; iloc++)
            {
                MATDTYPE* const localX_iloc = localX_->getSubMatrix(iloc);
                for (int icolor = 0; icolor < chromatic_number_; icolor++)
                    for (int jcolor = 0; jcolor < chromatic_number_; jcolor++)
                        localX_iloc[icolor + chromatic_number_ * jcolor]
                            = submatWork_->val(icolor, jcolor, iloc);
            }
        update_submatX_tm_.stop();
    }

    void updateSubMatLS() { gm_->updateSubMatLS(*submatLS_); }

    void updateSubMatT()
    {
        update_submatT_tm_.start();
        submatWork_->gather(*theta_);

        if (chromatic_number_ > 0)
            for (short iloc = 0; iloc < subdiv_; iloc++)
            {
                MATDTYPE* localT_iloc = localT_->getSubMatrix(iloc);
                for (int icolor = 0; icolor < chromatic_number_; icolor++)
                    for (int jcolor = 0; jcolor < chromatic_number_; jcolor++)
                        localT_iloc[icolor + chromatic_number_ * jcolor]
                            = submatWork_->val(icolor, jcolor, iloc);
            }
        update_submatT_tm_.stop();
    }

    void getLoewdinTransform(SquareLocalMatrices<MATDTYPE>& localP);

    void printTimers(ostream& os);

    void initializeGramMatrix(
        const SquareLocalMatrices<MATDTYPE>& ss, const int orbitals_index)
    {
        assert(gm_ != 0);

        init_gram_matrix_tm_.start();
        ss.fillDistMatrix(*work_, global_indexes_);

        gm_->setMatrix(*work_, orbitals_index);
        init_gram_matrix_tm_.stop();
    }
    void printGramMM(ofstream& tfile) { gm_->printMM(tfile); }
    void printHamiltonianMM(ofstream& tfile) { matH_->printMM(tfile); }
    int getDMMatrixIndex() const
    {
        assert(dm_ != 0);
        return dm_->getOrbitalsIndex();
    }
    void setDMuniform(const DISTMATDTYPE nel, const int orbitals_index)
    {
        dm_->setUniform(nel, orbitals_index);
    }
    int dim() const { return dim_; }

    bool withSpin() const { return with_spin_; }

    void computeInvS();

    virtual void computeInvB()
    {
        compute_invB_tm_.start();
        computeInvS();
        compute_invB_tm_.stop();
    }

    const dist_matrix::DistMatrix<DISTMATDTYPE>& dm() const;
    const dist_matrix::DistMatrix<DISTMATDTYPE>& kernel4dot() const;

    bool occupationsUptodate() const { return dm_->occupationsUptodate(); }

    const dist_matrix::DistMatrix<DISTMATDTYPE>& getLS() const
    {
        return gm_->getCholeskyL();
    }

    void updateLS() { gm_->updateLS(); }

    virtual void updateTheta()
    {
        // theta = invB * Hij
        update_theta_tm_.start();
        theta_->symm('l', 'l', 1., gm_->getInverse(), *matH_, 0.);

        update_theta_tm_.stop();
    }

    void updateThetaAndHB()
    {
        updateTheta();

        updateHB();
    }

    void assignH(const dist_matrix::DistMatrix<DISTMATDTYPE>& matH)
    {
        *matH_ = matH;
    }
    const dist_matrix::DistMatrix<DISTMATDTYPE>& getH() { return *matH_; }

    void setHB2H()
    {
        // if( onpe0 )
        //    (*MPIdata::sout)<<"ProjectedMatrices::setHB2H()..."<<endl;
        *matHB_ = *matH_;
    }

    virtual void updateHB() { setHB2H(); }

    void setGram2Id(const int orbitals_index) { gm_->set2Id(orbitals_index); }

    void setGramMatrix(const dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
        const int orbitals_index)
    {
        assert(gm_ != 0);
        gm_->setMatrix(mat, orbitals_index);
        gm_->computeInverse();
    }

    const dist_matrix::DistMatrix<DISTMATDTYPE>& getGramMatrix()
    {
        return gm_->getMatrix();
    }

    void getOccupations(std::vector<DISTMATDTYPE>& occ) const;
    void setOccupations(const std::vector<DISTMATDTYPE>& occ);

    void stripDM();
    void dressupDM();

    void applyInvS(SquareLocalMatrices<MATDTYPE>& mat);
    void applyInvS(dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
    {
        assert(gm_ != 0);
        gm_->applyInv(mat);
    }

    void setDMto2InvS();
    void buildDM(const dist_matrix::DistMatrix<DISTMATDTYPE>& z,
        const int orbitals_index);
    void buildDM(const dist_matrix::DistMatrix<DISTMATDTYPE>& z,
        const std::vector<DISTMATDTYPE>&, const int orbitals_index);
    void buildDM(const std::vector<DISTMATDTYPE>&, const int orbitals_index);

    double getEigSum();
    double getExpectation(const dist_matrix::DistMatrix<DISTMATDTYPE>& A);
    double getExpectationH();

    void solveGenEigenProblem(dist_matrix::DistMatrix<DISTMATDTYPE>& zz,
        std::vector<DISTMATDTYPE>& val, char job = 'v');
    void computeOccupationsFromDM();

    virtual void rotateAll(
        const dist_matrix::DistMatrix<DISTMATDTYPE>& rotation_matrix,
        const bool flag_eigen);

    double getNel() const;

    double getLowestEigenvalue() const { return eigenvalues_[0]; }

    void printS(ostream& os) const { gm_->print(os); }

    virtual void printMatrices(ostream& os) const
    {
        printS(os);
        printH(os);
        printTheta(os);
        printHB(os);
    }

    void printDM(ostream& os) const;
    void printOccupations(ostream& os) const;

    double computeCond() { return gm_->computeCond(); }

    double computeEntropy(const double kbt);
    double computeEntropy();
    double checkCond(const double tol, const bool flag = true);
    int writeDM_hdf5(HDFrestart& h5f_file);
    int read_dm_hdf5(hid_t file_id);
    void printEigenvalues(ostream& os) const;
    void setAuxilliaryEnergiesFromOccupations();
    void updateDM(const int iterative_index);
    void updateDMwithEigenstates(const int iterative_index);
    void updateDMwithEigenstatesAndRotate(
        const int iterative_index, dist_matrix::DistMatrix<DISTMATDTYPE>& zz);
    double computeChemicalPotentialAndOccupations(
        const double width, const int nel, const int max_numst)
    {
        return computeChemicalPotentialAndOccupations(
            aux_energies_, width, nel, max_numst);
    }

    dist_matrix::DistMatrix<DISTMATDTYPE> getDistMatrixFromLocalMatrices(
        const LocalMatrices<MATDTYPE>& ss)
    {
        dist_matrix::DistMatrix<DISTMATDTYPE> tmp("tmp", dim_, dim_);

        ss.fillDistMatrix(tmp, global_indexes_);

        return tmp;
    }

    void saveDM()
    {
        if (onpe0) cout << "ProjectedMatrices::saveDM()" << endl;
        if (!mat_X_old_)
        {
            mat_X_old_
                = new dist_matrix::DistMatrix<DISTMATDTYPE>(dm_->getMatrix());
            mat_L_old_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                gm_->getCholeskyL());
        }
        else
        {
            *mat_X_old_ = dm_->getMatrix();
            *mat_L_old_ = gm_->getCholeskyL();
        }
    }

    void resetDM()
    {
        dm_->setMatrix(*mat_X_old_, 0);
        dm_->stripS(*mat_L_old_, 0);
    }

    void updateDMwithRelax(const double mix, const int itindex)
    {
        // cout<<"ProjectedMatrices::updateDMwithRelax()..."<<endl;
        assert(mat_X_old_);

        dist_matrix::DistMatrix<DISTMATDTYPE> tmp(dm_->getMatrix());
        tmp.scal(mix);

        double alpha = 1. - mix;
        tmp.axpy(alpha, *mat_X_old_);
        dm_->setMatrix(tmp, itindex);
    }

    void getReplicatedDM(DISTMATDTYPE* replicated_DM_matrix)
    {
        const dist_matrix::DistMatrix<DISTMATDTYPE>& dm(dm_->getMatrix());

        dm.matgather(replicated_DM_matrix, dim_);
    }

    double getLinDependent2states(int& st1, int& st2, const bool flag) const
    {
        return gm_->getLinDependent2states(st1, st2);
    }

    static void registerRemoteTasksDistMatrix(
        dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>*
            remote_tasks_DistMatrix)
    {
        assert(remote_tasks_DistMatrix != 0);
        remote_tasks_DistMatrix_ = remote_tasks_DistMatrix;
    }
    virtual void initializeMatB(const SquareLocalMatrices<MATDTYPE>& ss)
    {
        (void)ss;
        cerr << "ERROR: ProjectedMatrices::initializeMatB() should never be "
                "called!!!"
             << endl;
    }
#if EXTRAPOLATE_H
    void initExtrapolationH() { *new_h_ = *matHB_; }
    void extrapolateHorder2(dist_matrix::DistMatrix<DISTMATDTYPE> matQ,
        dist_matrix::DistMatrix<DISTMATDTYPE>& yyt, ostream& os)
    {
        rotateSym(*h_minus1_, matQ, *work_);

        h_minus1_->axpy(-1., *new_h_);
        work_->gemm('n', 'n', 1., *h_minus1_, yyt, 0.);
        h_minus1_->gemm('t', 'n', 1., yyt, *work_, 0.);
        if (onpe0) os << "delta H..." << endl;
        h_minus1_->print(os, 0, 0, 5, 5);
        new_h_->axpy(-1., *h_minus1_);
    }
    void updateHminus1tmp(dist_matrix::DistMatrix<DISTMATDTYPE> matQ,
        dist_matrix::DistMatrix<DISTMATDTYPE>& yyt, ostream& os)
    {
        if (h_minus1_ != 0)
        {
            rotateSym(*h_minus1_, matQ, *work_);
            *tmp_h_minus1_ = *h_minus1_;
            tmp_h_minus1_->axpy(-1., *new_h_);
            work_->gemm('n', 'n', 1., *tmp_h_minus1_, yyt, 0.);
            tmp_h_minus1_->gemm('t', 'n', 1., yyt, *work_, 0.);
            if (onpe0) os << "delta H..." << endl;
            tmp_h_minus1_->print(os, 0, 0, 5, 5);
        }
    }
    void updateHminus2(dist_matrix::DistMatrix<DISTMATDTYPE> matQ,
        dist_matrix::DistMatrix<DISTMATDTYPE>& yyt, ostream& os)
    {
        if (h_minus2_ != 0)
        {
            rotateSym(*h_minus2_, matQ, *work_);
            h_minus2_->axpy(-1., *h_minus1_);
            work_->gemm('n', 'n', 1., *h_minus2_, yyt, 0.);
            h_minus2_->gemm('t', 'n', 1., yyt, *work_, 0.);
            if (onpe0) os << "delta H..." << endl;
            h_minus2_->print(os, 0, 0, 5, 5);
        }
    }
    void extrapolateHorder3()
    {
        if (h_minus2_ != 0)
        {
            new_h_->axpy(-2., *tmp_h_minus1_);
            new_h_->axpy(1., *h_minus2_);
        }
        else
        {
            if (h_minus1_ != 0)
            {
                new_h_->axpy(-1., *tmp_h_minus1_);
            }
        }
    }

    void saveH() { *matHB_ = *new_h_; }
    void updateHminus1()
    {
        if (h_minus1_ == 0)
        {
            h_minus1_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                "H_minus1", dim_, dim_);
        }
        *h_minus1_ = *matHB_;
    }
    void updateHminus2()
    {
        if (h_minus1_ == 0)
        {
            h_minus1_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                "H_minus1", dim_, dim_);
        }
        else
        {
            if (h_minus2_ == 0)
                h_minus2_ = new dist_matrix::DistMatrix<DISTMATDTYPE>(
                    "H_minus2", dim_, dim_);
            *h_minus2_ = *h_minus1_;
        }
        *h_minus1_ = *matHB_;
    }
#endif
#if 0
    void rotateBackDM();
#endif

    void resetDotProductMatrices();
    double dotProductWithInvS(const SquareLocalMatrices<MATDTYPE>& ss);
    double dotProductWithDM(const SquareLocalMatrices<MATDTYPE>& ss);
    double dotProductSimple(const SquareLocalMatrices<MATDTYPE>& ss);
    double computeTraceInvSmultMat(const SquareLocalMatrices<MATDTYPE>& mat);
    double computeTraceInvSmultMat(
        const dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
    {
        dist_matrix::DistMatrix<DISTMATDTYPE> pmatrix(mat);
        gm_->applyInv(pmatrix);

        return pmatrix.trace();
    }
    double computeTraceInvSmultMatMultTheta(
        const dist_matrix::DistMatrix<DISTMATDTYPE>& mat);
    double computeTraceMatMultTheta(
        const dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
    {
        assert(theta_ != 0);
        dist_matrix::DistMatrix<DISTMATDTYPE> pmatrix("pmatrix", dim_, dim_);
        pmatrix.gemm('n', 'n', 1.0, mat, *theta_, 0.);
        return pmatrix.trace();
    }
    void computeMatMultTheta(const dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
        dist_matrix::DistMatrix<DISTMATDTYPE>& pmat)
    {
        assert(theta_ != 0);
        pmat.gemm('n', 'n', 1.0, mat, *theta_, 0.);
    }
    const dist_matrix::DistMatrix<DISTMATDTYPE>& getMatHB() { return *matHB_; }
    void setDM(const dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
        const int orbitals_index)
    {
        dm_->setMatrix(mat, orbitals_index);
    }
    void computeGenEigenInterval(std::vector<double>& interval,
        const int maxits, const double padding = 0.01);
};

#endif
