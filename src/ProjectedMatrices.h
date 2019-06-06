// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
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
#include "DistMatrix2SquareLocalMatrices.h"
#include "LocalMatrices2DistMatrix.h"

#include "hdf5.h"

#include <iostream>


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
    std::unique_ptr<dist_matrix::DistMatrix<DISTMATDTYPE>> mat_X_old_;
    std::unique_ptr<dist_matrix::DistMatrix<DISTMATDTYPE>> mat_L_old_;

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

    ProjectedMatrices& operator=(const ProjectedMatrices& src);
    ProjectedMatrices(const ProjectedMatrices& pm);

    /*!
     * Matrices to be used to multiply orbitals on the right
     * (dimension equal to number of local colors)
     */
    std::unique_ptr<SquareLocalMatrices<MATDTYPE>> localX_; // density matrix
    std::unique_ptr<SquareLocalMatrices<MATDTYPE>> localT_; // theta=inv(S)*H_phi

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

    std::unique_ptr<dist_matrix::DistMatrix<DISTMATDTYPE>> theta_;
    std::unique_ptr<dist_matrix::DistMatrix<DISTMATDTYPE>> matHB_;
    std::unique_ptr<dist_matrix::DistMatrix<DISTMATDTYPE>> matH_;

    std::unique_ptr<DensityMatrix> dm_;
    std::unique_ptr<GramMatrix> gm_;

    // work matrix for tmp usage
    std::unique_ptr<dist_matrix::DistMatrix<DISTMATDTYPE>> work_;

    void printTheta(std::ostream& os) const
    {
        if (onpe0) os << " Matrix Theta" << std::endl;
        theta_->print(os, 0, 0, NPRINT_ROWS_AND_COLS, NPRINT_ROWS_AND_COLS);
    }

    void printHB(std::ostream& os) const
    {
        if (onpe0) os << " Matrix HB" << std::endl;
        matHB_->print(os, 0, 0, NPRINT_ROWS_AND_COLS, NPRINT_ROWS_AND_COLS);
    }

    void printH(std::ostream& os) const
    {
        if (onpe0) os << " Matrix H" << std::endl;
        matH_->print(os, 0, 0, NPRINT_ROWS_AND_COLS, NPRINT_ROWS_AND_COLS);
    }

public:
    ProjectedMatrices(const int, const bool with_spin);
    virtual ~ProjectedMatrices();

    virtual void setup(const double kbt, const int nel,
        const std::vector<std::vector<int>>& global_indexes);

    void addMatrixElementSparseH(const int st1, const int st2, const double val)
    {
        sH_->push_back(st1, st2, val);
    }

    // fill SparseDistMatrix sH_ with values in slH
    void addMatrixElementsSparseH(const SquareLocalMatrices<MATDTYPE>& slH)
    {
        Control& ct = *(Control::instance());

       LocalMatrices2DistMatrix* sl2dm =
            LocalMatrices2DistMatrix::instance();

        sl2dm->convert(slH, *sH_, ct.numst);
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

    const dist_matrix::SubMatrices<DISTMATDTYPE>& getSubMatLS(
        MPI_Comm comm, const std::vector<std::vector<int>>& global_indexes) const
    {
        return gm_->getSubMatLS(comm, global_indexes);
    }

    void updateSubMatX() { updateSubMatX(dm_->getMatrix()); }

    void updateSubMatX(const dist_matrix::DistMatrix<DISTMATDTYPE>& dm)
    {
        DistMatrix2SquareLocalMatrices* dm2sl =
            DistMatrix2SquareLocalMatrices::instance();
        dm2sl->convert(dm, *localX_);
    }

    void updateSubMatT()
    {
        DistMatrix2SquareLocalMatrices* dm2sl =
            DistMatrix2SquareLocalMatrices::instance();
        dm2sl->convert(*theta_, *localT_);
    }

    void getLoewdinTransform(SquareLocalMatrices<MATDTYPE>& localP);

    void printTimers(std::ostream& os);

    void initializeGramMatrix(
        const SquareLocalMatrices<MATDTYPE>& ss, const int orbitals_index)
    {
        assert(gm_);

        init_gram_matrix_tm_.start();

        LocalMatrices2DistMatrix* sl2dm =
            LocalMatrices2DistMatrix::instance();

        sl2dm->convert(ss, *work_, dim_);

        gm_->setMatrix(*work_, orbitals_index);
        init_gram_matrix_tm_.stop();
    }
    void printGramMM(std::ofstream& tfile) { gm_->printMM(tfile); }
    void printHamiltonianMM(std::ofstream& tfile) { matH_->printMM(tfile); }
    int getDMMatrixIndex() const
    {
        assert(dm_);
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

    void printS(std::ostream& os) const { gm_->print(os); }

    virtual void printMatrices(std::ostream& os) const
    {
        printS(os);
        printH(os);
        printTheta(os);
        printHB(os);
    }

    void printDM(std::ostream& os) const;
    void printOccupations(std::ostream& os) const;

    double computeCond() { return gm_->computeCond(); }

    double computeEntropy(const double kbt);
    double computeEntropy();
    double checkCond(const double tol, const bool flag = true);
    int writeDM_hdf5(HDFrestart& h5f_file);
    int read_dm_hdf5(hid_t file_id);
    void printEigenvalues(std::ostream& os) const;
    void updateDM(const int iterative_index);
    void updateDMwithEigenstates(const int iterative_index);
    void updateDMwithEigenstatesAndRotate(
        const int iterative_index, dist_matrix::DistMatrix<DISTMATDTYPE>& zz);
    double computeChemicalPotentialAndOccupations(
        const double width, const int nel, const int max_numst)
    {
        return computeChemicalPotentialAndOccupations(
            eigenvalues_, width, nel, max_numst);
    }

    dist_matrix::DistMatrix<DISTMATDTYPE> getDistMatrixFromLocalMatrices(
        const LocalMatrices<MATDTYPE>& ss)
    {
        dist_matrix::DistMatrix<DISTMATDTYPE> tmp("tmp", dim_, dim_);

        LocalMatrices2DistMatrix* sl2dm =
            LocalMatrices2DistMatrix::instance();

        sl2dm->convert(ss, tmp, dim_);

        return tmp;
    }

    void saveDM()
    {
        if (onpe0) std::cout << "ProjectedMatrices::saveDM()" << std::endl;
        if (!mat_X_old_)
        {
            mat_X_old_.reset(
                new dist_matrix::DistMatrix<DISTMATDTYPE>(dm_->getMatrix()));
            mat_L_old_.reset( new dist_matrix::DistMatrix<DISTMATDTYPE>(
                gm_->getCholeskyL()) );
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
        std::cerr << "ERROR: ProjectedMatrices::initializeMatB() should never be "
                "called!!!"
             << std::endl;
    }

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
    dist_matrix::DistMatrix<DISTMATDTYPE>& getMatHB() { return *matHB_; }
    void setDM(const dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
        const int orbitals_index)
    {
        dm_->setMatrix(mat, orbitals_index);
    }
    void computeGenEigenInterval(std::vector<double>& interval,
        const int maxits, const double padding = 0.01);
    DensityMatrix& getDM() { return *dm_; }
};

#endif
