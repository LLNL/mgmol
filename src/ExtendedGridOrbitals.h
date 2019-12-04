// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_EXTENDEDGRIDORBITALS_H
#define MGMOL_EXTENDEDGRIDORBITALS_H

#include "BlockVector.h"
#include "DistMatrix.h"
#include "GridFunc.h"
#include "HDFrestart.h"
#include "Lap.h"
#include "MPIdata.h"
#include "Mesh.h"
#include "Orbitals.h"
#include "SinCosOps.h"
#include "global.h"

#include "hdf5.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

template <class T>
class LocalMatrices;
template <class T>
class SquareLocalMatrices;
class Potentials;
class ProjectedMatrices;
class ProjectedMatricesInterface;
class LocalizationRegions;
class ExtendedGridOrbitals;
class MasksSet;
class ClusterOrbitals;

typedef double (ExtendedGridOrbitals::*ExtendedGridOrbitalsPtrFunc)(
    const ExtendedGridOrbitals&);

class ExtendedGridOrbitals : public Orbitals
{
private:
    const std::string name_;

    ////////////////////////////////////////////////////////
    // common data shared by all instances of class
    ////////////////////////////////////////////////////////
    static Timer matB_tm_;
    static Timer invBmat_tm_;
    static Timer overlap_tm_;
    static Timer dot_product_tm_;
    static Timer addDot_tm_;
    static Timer prod_matrix_tm_;
    static Timer assign_tm_;
    static Timer normalize_tm_;
    static Timer axpy_tm_;

    static short bc_[3];

    static int lda_; // leading dimension for storage
    static int numpt_;
    static int loc_numpt_;

    // static double (ExtendedGridOrbitals::*dotProduct_)(const
    // ExtendedGridOrbitals&);
    static ExtendedGridOrbitalsPtrFunc dotProduct_;

    static int data_wghosts_index_;

    static int numst_;

    ////////////////////////////////////////////////////////
    // common data shared by copies of object (and copied whan a copy is made)
    ////////////////////////////////////////////////////////

    // pointers to objects owned outside class
    ProjectedMatricesInterface* proj_matrices_;

    ////////////////////////////////////////////////////////
    // instance specific data
    ////////////////////////////////////////////////////////
    BlockVector<ORBDTYPE, MemorySpace::Host> block_vector_;

    ////////////////////////////////////////////////////////
    //
    // private functions
    //
    void projectOut(ORBDTYPE* const, const int, const double scale = 1.);

    void multiply_by_matrix(
        const DISTMATDTYPE* const, ORBDTYPE*, const int) const;
    void multiply_by_matrix(const dist_matrix::DistMatrix<DISTMATDTYPE>& matrix,
        ORBDTYPE* const product, const int ldp);
    void scal(const int i, const double alpha) { block_vector_.scal(i, alpha); }
    virtual void assign(const int i, const ORBDTYPE* const v, const int n = 1)
    {
        block_vector_.assign(i, v, n);
    }
    ExtendedGridOrbitals& operator=(const ExtendedGridOrbitals& orbitals);
    ExtendedGridOrbitals();

    void computeMatB(const ExtendedGridOrbitals&, const pb::Lap<ORBDTYPE>&);

    double dotProductDiagonal(const ExtendedGridOrbitals& orbitals);
    double dotProductWithDM(const ExtendedGridOrbitals& orbitals);
    double dotProductWithInvS(const ExtendedGridOrbitals& orbitals);
    double dotProductSimple(const ExtendedGridOrbitals& orbitals);

    void computeLocalProduct(const ORBDTYPE* const, const int,
        LocalMatrices<MATDTYPE>&, const bool transpose = false);

    void computeGlobalIndexes();
    void computeInvNorms2(std::vector<std::vector<double>>& inv_norms2) const;
    void computeDiagonalGram(VariableSizeMatrix<sparserow>& diagS) const;

    void initFourier();
    void initRand();
    dist_matrix::DistMatrix<DISTMATDTYPE> product(const ORBDTYPE* const,
        const int, const int, const bool transpose = false);

    ORBDTYPE* psi(const int i) const { return block_vector_.vect(i); }

    void app_mask(const int, ORBDTYPE*, const short) const {};
    void multiplyByMatrix(const SquareLocalMatrices<MATDTYPE>& matrix,
        ORBDTYPE* product, const int ldp) const;
    void setup();

protected:
    const pb::Grid& grid_;

    static short subdivx_;

    // indexes corresponding to valid function in each subdomain
    static std::vector<std::vector<int>> overlapping_gids_;

public:
    friend class SinCosOps<ExtendedGridOrbitals>;

    double norm() const;

    ProjectedMatricesInterface* proj_matrices() const { return proj_matrices_; }

    const std::vector<int>& getAllOverlappingGids() const
    {
        return overlapping_gids_[0];
    }
    const std::vector<int>& getLocalGids() const
    {
        return overlapping_gids_[0];
    }

    ExtendedGridOrbitals(std::string name, const pb::Grid& my_grid,
        const short subdivx, const int numst, const short bc[3],
        ProjectedMatricesInterface*, LocalizationRegions*, MasksSet* masks,
        MasksSet* corrmasks, ClusterOrbitals* local_cluster,
        const bool setup_flag = true);

    ExtendedGridOrbitals(const std::string& name, const ExtendedGridOrbitals& A,
        const bool copy_data = true);
    ExtendedGridOrbitals(const std::string& name, const ExtendedGridOrbitals& A,
        ProjectedMatricesInterface* proj_matrices, const bool copy_data = true);

    virtual ~ExtendedGridOrbitals();

    static void printTimers(std::ostream& os);

    void resetDotProductMatrices();
    void init2zero();

    void reset(MasksSet* masks, MasksSet* corrmasks, LocalizationRegions* lrs);

    virtual void assign(const ExtendedGridOrbitals& orbitals);
    void copyDataFrom(const ExtendedGridOrbitals& src);

    ProjectedMatricesInterface* getProjMatrices() { return proj_matrices_; }

    const ProjectedMatricesInterface* projMatrices() const
    {
        return proj_matrices_;
    }

    int numst(void) const { return numst_; }
    int getLda() const { return lda_; }
    int getLocNumpt() const { return loc_numpt_; }

    bool isCompatibleWith(const ExtendedGridOrbitals&) const { return true; }

    void resetDataWithGhostsIndex() { data_wghosts_index_ = 0; }

    void setDataWithGhosts(const bool force = false)
    {
        if (data_wghosts_index_ == getIterativeIndex() && !force) return;

        block_vector_.setDataWithGhosts();

        // if( onpe0 )
        //    (*MPIdata::sout)<<"setDataWithGhosts with iterative index
        //    "<<getIterativeIndex()<<endl;
        data_wghosts_index_ = getIterativeIndex();
    }

    template <typename T>
    void setDataWithGhosts(pb::GridFuncVector<T>* data_wghosts)
    {
        assert(data_wghosts != 0);

        block_vector_.setDataWithGhosts(data_wghosts);
    }
    pb::GridFunc<ORBDTYPE>& getFuncWithGhosts(const int i)
    {
        //(*MPIdata::sout)<<" data_wghosts_index_="<<data_wghosts_index_
        //    <<" getIterativeIndex()   ="<<getIterativeIndex()<<endl;
        if (data_wghosts_index_ != getIterativeIndex())
        {
            setDataWithGhosts();
        }
        assert(data_wghosts_index_ == getIterativeIndex());
        //(*MPIdata::sout)<<"getFuncWithGhosts with index
        //"<<getIterativeIndex()<<endl;
        return block_vector_.getVectorWithGhosts(i);
    }

    pb::GridFuncVector<ORBDTYPE>* getPtDataWGhosts()
    {
        return block_vector_.getPtDataWGhosts();
    }

    void trade_boundaries()
    {
        static int last_index_traded = -1;

        if (data_wghosts_index_ != last_index_traded)
        {
#ifdef PRINT_OPERATIONS
            if (onpe0)
                (*MPIdata::sout)
                    << "ExtendedGridOrbitals::trade_boundaries()" << std::endl;
#endif
            block_vector_.trade_boundaries();
            last_index_traded = data_wghosts_index_;
        }
    }

    void set_storage(ORBDTYPE* new_storage)
    {
        assert(new_storage != 0);
        block_vector_.setStorage(new_storage);
    }
    ORBDTYPE* getPsi(const int i, const short iloc = 0) const
    {
        assert(iloc < subdivx_);
        return block_vector_.vect(i) + iloc * loc_numpt_;
    }
    template <typename T>
    void setPsi(const pb::GridFunc<T>& gf_work, const int ist)
    {
        block_vector_.assignComponent(gf_work, ist);
    }
    template <typename T>
    void setPsi(const pb::GridFuncVector<T>& gf_work)
    {
        block_vector_.assign(gf_work);
    }
    int chromatic_number(void) const
    {
        assert(numst_ < 10000);
        return numst_;
    }
    short subdivx(void) const { return subdivx_; }
    void printChromaticNumber(std::ostream& os) const
    {
        if (onpe0) os << " Max. chromatic_number: " << numst_ << std::endl;
    }
    void printNumst(std::ostream& os) const
    {
        if (onpe0) os << " Number of states   = " << numst_ << std::endl;
    }
    void computeBAndInvB(const pb::Lap<ORBDTYPE>& LapOper);

    void computeGram(const int verbosity = 0);
    void computeGramAndInvS(const int verbosity = 0);
    void computeGram(dist_matrix::DistMatrix<DISTMATDTYPE>& gram_mat);
    void computeGram(const ExtendedGridOrbitals& orbitals,
        dist_matrix::DistMatrix<DISTMATDTYPE>& gram_mat);

    ORBDTYPE maxAbsValue() const { return block_vector_.maxAbsValue(); }

    /*!
     * use predefined (default) dot product type
     */
    double dotProduct(const ExtendedGridOrbitals& orbitals);
    /*!
     * use different dot product type
     */
    double dotProduct(const ExtendedGridOrbitals&, const short dot_type);

    static void setDotProduct(const short dot_type);
    void computeDiagonalElementsDotProduct(const ExtendedGridOrbitals& orbitals,
        std::vector<DISTMATDTYPE>& ss) const;

    dist_matrix::DistMatrix<DISTMATDTYPE> product(
        const ExtendedGridOrbitals&, const bool transpose = false);
    void computeLocalProduct(const ExtendedGridOrbitals&,
        LocalMatrices<MATDTYPE>&, const bool transpose = false);
    void getLocalOverlap(SquareLocalMatrices<MATDTYPE>&);
    void getLocalOverlap(
        const ExtendedGridOrbitals& orbitals, SquareLocalMatrices<MATDTYPE>&);

    void addDotWithNcol2Matrix(
        ExtendedGridOrbitals&, dist_matrix::DistMatrix<double>&) const;

    void scal(const double alpha)
    {
        block_vector_.scal(alpha);
        incrementIterativeIndex();
    }
    void projectOut(ExtendedGridOrbitals&, const double scale = 1.);

    void normalize();
    void orthonormalize2states(const int st1, const int st2);
    void orthonormalizeLoewdin(const bool overlap_uptodate = false,
        SquareLocalMatrices<MATDTYPE>* matrixTransform     = nullptr,
        const bool update_matrices                         = true);

    ExtendedGridOrbitals& operator-=(const ExtendedGridOrbitals& orbitals)
    {
        block_vector_ -= orbitals.block_vector_;
        return *this;
    }

    void initGauss(const double, const LocalizationRegions&);
    virtual void axpy(const double alpha, const ExtendedGridOrbitals&);

    void app_mask(const int, pb::GridFunc<ORBDTYPE>&, const short) const {};

    void applyMask(const bool = false){};
    void applyCorrMask(const bool = false){};

    void multiplyByMatrix(const SquareLocalMatrices<MATDTYPE>& matrix);
    void multiplyByMatrix(const SquareLocalMatrices<MATDTYPE>& matrix,
        ExtendedGridOrbitals& product) const;
    void multiply_by_matrix(
        const DISTMATDTYPE* const matrix, ExtendedGridOrbitals& product) const;
    void multiply_by_matrix(const dist_matrix::DistMatrix<DISTMATDTYPE>&);
    void multiplyByMatrix2states(const int st1, const int st2,
        const double* mat, ExtendedGridOrbitals& product);

    int write_hdf5(HDFrestart& h5f_file, const std::string& name = "Function");
    int write_func_hdf5(HDFrestart&, const std::string& name = "Function");
    int read_hdf5(HDFrestart& h5f_file);
    int read_func_hdf5(HDFrestart&, const std::string& name = "Function");

    void initWF(const LocalizationRegions& lrs);
    void checkCond(const double tol, const bool flag_stop);
    double normState(const int st) const;
    const std::vector<std::vector<int>>& getOverlappingGids() const
    {
        assert(overlapping_gids_.size() > 0);
        return overlapping_gids_;
    }
    int getGlobalIndex(const short iloc, const short color) const
    {
        assert(overlapping_gids_.size() > 0);
        assert(iloc < static_cast<int>(overlapping_gids_.size()));
        assert(color < static_cast<int>(overlapping_gids_[iloc].size()));
        return overlapping_gids_[iloc][color];
    }
    int getColor(const int gid) const { return gid; }
    double getMaxR() const
    {
        Mesh* mymesh           = Mesh::instance();
        const pb::Grid& mygrid = mymesh->grid();
        return mygrid.maxDomainSize();
    }
};

#endif
