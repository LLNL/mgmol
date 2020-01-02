// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_PROJECTED_MATRICES_INTERFACE_H
#define MGMOL_PROJECTED_MATRICES_INTERFACE_H

#include <cassert>
#include <iostream>
#include <vector>

/* Currently only needed functions in ProjectedMatrices.h class. Must be
 * cleaned-up later */
#include "DistMatrix.h"
#include "MGmol_MPI.h"
#include "SquareLocalMatrices.h"
#include "SquareSubMatrix.h"
#include "VariableSizeMatrix.h"
#include "entropy.h"
#include "fermi.h"
#include "hdf5.h"
#include "tools.h"

typedef DISTMATDTYPE PROJMATDTYPE;

class HDFrestart;

/* tolerance for matrix elements */
const double tol_matrix_elements = 1.e-14;

class ProjectedMatricesInterface
{
private:
    int iterative_index_h_;

    int index(const int phi_index, const int pot_index) const
    {
        return (phi_index % 100) + 100 * pot_index;
    }

protected:
    short subdiv_;
    short chromatic_number_;
    double width_;
    double mu_;
    int nel_;

public:
    ProjectedMatricesInterface()
        : iterative_index_h_(-1),
          subdiv_(-1),
          chromatic_number_(-1),
          width_(0.),
          nel_(0.){};

    virtual ~ProjectedMatricesInterface(){};

    virtual void setup(const double kbt, const int nel,
        const std::vector<std::vector<int>>& global_indexes)
        = 0;

    // initial setup function
    void setupBase(const double kbt, const int num_el, const int subdiv,
        const int chromatic_number)
    {
        width_            = kbt;
        nel_              = num_el;
        subdiv_           = subdiv;
        chromatic_number_ = chromatic_number;
    }

    void setHiterativeIndex(const int phi_index, const int pot_index)
    {
        assert(phi_index < 100000);
        iterative_index_h_ = index(phi_index, pot_index);
    }
    bool isHupToDate(const int phi_index, const int pot_index)
    {
        assert(phi_index < 100000);
        return (iterative_index_h_ == index(phi_index, pot_index));
    }
    virtual void printTimers(std::ostream& os) = 0;
    // initialize Gram matrix with SquareLocalMatrices (input)
    virtual void initializeGramMatrix(
        const SquareLocalMatrices<MATDTYPE>& ss, const int orbitals_index)
        = 0;

    virtual void setLocalMatrixElementsHl(
        const SquareLocalMatrices<MATDTYPE>& slH)
        = 0;
    virtual void setLocalMatrixElementsHnl(const SquareSubMatrix<MATDTYPE>& slH)
        = 0;

    virtual void clearSparseH() = 0;
    virtual void consolidateH() = 0;

    // return density matrix (inverse Gram if no unoccupied states)
    virtual SquareLocalMatrices<MATDTYPE>& getLocalX() const = 0;

    // return S**(-1)*H (or B**(-1)*H with Mehrstellen)
    virtual SquareLocalMatrices<MATDTYPE>& getLocalT() const = 0;

    // returns local Gram Matrix
    //    virtual SquareLocalmatrices& getLocalS()const=0;

    // return Tr(X*H) or Tr(S**(-1)*H) if occupied states only
    virtual double getExpectationH() = 0;

    virtual void updateSubMatX() = 0;

    virtual void updateSubMatT() = 0;

    // return sum_i ( ddiag[i]*S**(-1)[i][i] )
    virtual double getTraceDiagProductWithInvS(std::vector<PROJMATDTYPE>& ddiag)
        = 0;

    virtual double computeEntropy() = 0;
    virtual double getEigSum()      = 0;

    virtual void updateTheta()                     = 0;
    virtual void computeInvB()                     = 0;
    virtual void printGramMM(std::ofstream& tfile) = 0;
    virtual void setDMuniform(const PROJMATDTYPE nel, const int orbitals_index)
        = 0;

    virtual double dotProductWithInvS(const SquareLocalMatrices<MATDTYPE>& ss)
        = 0;
    virtual double dotProductWithDM(const SquareLocalMatrices<MATDTYPE>& ss)
        = 0;
    virtual double dotProductSimple(const SquareLocalMatrices<MATDTYPE>& ss)
        = 0;
    virtual void resetDotProductMatrices() = 0;

    virtual void stripDM()                                     = 0;
    virtual void dressupDM()                                   = 0;
    virtual void printS(std::ostream& os) const                = 0;
    virtual void computeInvS()                                 = 0;
    virtual void printMatrices(std::ostream& os) const         = 0;
    virtual void applyInvS(SquareLocalMatrices<MATDTYPE>& mat) = 0;
    virtual double getLinDependent2states(
        int& st1, int& st2, const bool) const                          = 0;
    virtual double checkCond(const double tol, const bool flag = true) = 0;
    virtual double getNel() const                                      = 0;
    virtual void updateThetaAndHB()                                    = 0;
    virtual void printDM(std::ostream& os) const                       = 0;

    //////////////////////////////////////////////////////////////////////
    /* Default implementation - returns error message */
    virtual void printOccupations(std::ostream& os) const
    {
        (void)os;

        exitWithErrorMessage("printOccupations");
    }
    virtual int getDMMatrixIndex() const
    {
        exitWithErrorMessage("getDMMatrixIndex");

        return 0;
    }
    virtual double computeCond() = 0;
    virtual void getOccupations(std::vector<PROJMATDTYPE>& occ) const
    {
        (void)occ;

        exitWithErrorMessage("getOccupations");
    }
    virtual void setHB2H() { exitWithErrorMessage("setHB2H"); }
    virtual void setAuxilliaryEnergiesFromOccupations()
    {
        exitWithErrorMessage("setAuxilliaryEnergiesFromOccupations");
    }
    virtual const dist_matrix::DistMatrix<DISTMATDTYPE>& dm() const
    {
        exitWithErrorMessage("dm");

        dist_matrix::DistMatrix<DISTMATDTYPE>* tmp
            = new dist_matrix::DistMatrix<DISTMATDTYPE>("tmp");

        return (*tmp);
    }
    virtual void saveDM() { exitWithErrorMessage("saveDM"); }
    virtual void resetDM() { exitWithErrorMessage("resetDM"); }
    virtual void updateDMwithRelax(const double mix, const int itindex)
    {
        (void)mix;
        (void)itindex;

        exitWithErrorMessage("updateDMwithRelax");
    }
    virtual int read_dm_hdf5(hid_t file_id)
    {
        (void)file_id;

        exitWithErrorMessage("read_dm_hdf5");

        return 0;
    }
    virtual int writeDM_hdf5(HDFrestart& h5f_file)
    {
        (void)h5f_file;

        exitWithErrorMessage("writeDM_hdf5");

        return 0;
    }
    virtual void updateDM(const int iterative_index)
    {
        (void)iterative_index;

        exitWithErrorMessage("updateDM");
    }

    virtual void setDMto2InvS() { exitWithErrorMessage("setDMto2InvS"); }
    virtual void initializeMatB(const SquareLocalMatrices<MATDTYPE>& ss)
    {
        (void)ss;

        exitWithErrorMessage("initializeMatB");
    }
    virtual double getLowestEigenvalue() const = 0;
};

#endif
