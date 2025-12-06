// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef BLOCKVECTOR_H
#define BLOCKVECTOR_H

#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "memory_space.h"
#include "mputils.h"

#include "GridFuncVector.h"

#include <vector>

#define PAD 0

// When using NEWSTORAGE memory allocation is done for the class and not for
// specific instances thus the use of static variables and functions this is
// done to avoid not finding blocks large enough later during run due to memory
// fragmentation
template <typename ScalarType, typename MemorySpaceType>
class BlockVector
{
    static Timer set_data_tm_;
    static Timer trade_data_tm_;
    static Timer assign_tm_;
    static Timer scal_tm_;
    static Timer opminus_tm_;
    static Timer copy_tm_;
    static Timer diagop_tm_;

    static short n_instances_;
    static short subdivx_;

    static pb::GridFuncVector<ScalarType, MemorySpaceType>* data_wghosts_;

    // data allocator
    static std::vector<ScalarType*> class_storage_;

    // tells which blocks are allocated (1) or not (0)
    // to a specific instantiated object
    static std::vector<short> allocated_;

    static short max_alloc_instances_;
    static int size_storage_;
    static int allocated_size_storage_;
    int size_storage_instance_;

    // safety coefficient to multiply dize of initial allocation
    // and have enough memory for later allocations
    static float overallocate_factor_;

    // index telling which slot of memory has been allocated
    // for that object
    short my_allocation_;

    ScalarType* storage_; // storage for data

    std::vector<ScalarType*> vect_; // pointers to storage of vectors

    static int numel_; // nb elements in each vector
    static int locnumel_;
    static int ld_; // leading dimension

    int ld_instance_;

    const pb::Grid& mygrid_;
    short bc_[3];

    void allocate_storage();
    void deallocate_storage();
    BlockVector(const BlockVector& bv);
    void setup(const BlockVector& bv);

    static void allocate1NewBlock();

public:
    using memory_space_type = MemorySpaceType;

    BlockVector(
        const pb::Grid& my_grid, const short subdivx, const short bc[3]);

    BlockVector(const BlockVector& bv, const bool copy_data);

    BlockVector& operator=(const BlockVector& bv);

    ~BlockVector();

    const pb::GridFuncVector<ScalarType, MemorySpaceType>& getDataWGhosts()
    {
        assert(data_wghosts_ != 0);
        return *data_wghosts_;
    }

    pb::GridFuncVector<ScalarType, MemorySpaceType>* getPtDataWGhosts()
    {
        return data_wghosts_;
    }

    void initialize(
        const std::vector<std::vector<int>>& gid, const bool skinny_stencil);

    void clear()
    {
        vect_.clear();

        deallocate_storage();
    }

    template <typename ScalarType2>
    void axpy(const ScalarType2 alpha, const BlockVector& bv)
    {
        assert(storage_ != nullptr);
        assert(bv.storage_ != nullptr);
        assert(storage_ != bv.storage_);

        LinearAlgebraUtils<MemorySpaceType>::MPaxpy(
            size_storage_, alpha, bv.storage_, storage_);
    }

    void axpy(const double alpha, const ScalarType* const vy)
    {
        LinearAlgebraUtils<MemorySpaceType>::MPaxpy(
            size_storage_, alpha, vy, storage_);
    }

    ScalarType* vect(const int i) const
    {
        assert(i < (int)vect_.size());
        assert(vect_[i] != 0);
        return vect_[i];
    }

    void applyDiagonalOp(const std::vector<double>& diag,
        BlockVector<ScalarType, MemorySpaceType>& dst) const;

    ScalarType maxAbsValue() const;

    template <typename ScalarType2>
    void setDataWithGhosts(
        pb::GridFuncVector<ScalarType2, MemorySpaceType>* data_wghosts);

    void setDataWithGhosts();

    void assign(const int color, const ScalarType* const src, const int n = 1)
    {
        assert((color + n - 1) < static_cast<int>(vect_.size()));
        int ione   = 1;
        int mysize = n * ld_;
        Tcopy(&mysize, src, &ione, vect_[color], &ione);
    }

    void assign(const ScalarType* const src)
    {
        MemorySpace::Memory<ScalarType, MemorySpaceType>::copy(
            src, size_storage_, storage_);
    }

    /*
     * assign functions for source data with ghost values
     */
    template <typename ScalarType2>
    void assign(const pb::GridFuncVector<ScalarType2, MemorySpaceType>& src);
    template <typename ScalarType2>
    void assignComponent(const pb::GridFunc<ScalarType2>& src, const int i);
    template <typename ScalarType2>
    void assignComponent(
        const pb::GridFuncVector<ScalarType2, MemorySpaceType>& src,
        const int i);

    void assignLocal(
        const int color, const short iloc, const ScalarType* const src)
    {
        assert(color >= 0);
        assert(color < static_cast<int>(vect_.size()));
        assert(iloc < subdivx_);
        MemorySpace::Memory<ScalarType, MemorySpaceType>::copy(
            src, locnumel_, vect_[color] + iloc * locnumel_);
    }

    void setToDataWithGhosts() { assign(*data_wghosts_); }

    void copyDataFrom(const BlockVector& src)
    {
        copy_tm_.start();

        assert(src.size_storage_ == size_storage_);
        assert(storage_ != nullptr);
        assert(src.storage_ != nullptr);

        MemorySpace::Memory<ScalarType, MemorySpaceType>::copy(
            src.storage_, size_storage_, storage_);

        copy_tm_.stop();
    }

    pb::GridFunc<ScalarType>& getVectorWithGhosts(const int i)
    {
        assert(data_wghosts_ != 0);
        assert(i < static_cast<int>(vect_.size()));

        return data_wghosts_->getGridFunc(i);
    }

    void setStorage(ScalarType* new_storage)
    {
        assert(new_storage != 0 || vect_.size() == 0);
        if (std::is_same<MemorySpaceType, MemorySpace::Host>::value)
            MemorySpace::assert_is_host_ptr(new_storage);
        else
            MemorySpace::assert_is_dev_ptr(new_storage);

        storage_ = new_storage;
        for (unsigned int i = 0; i < vect_.size(); i++)
        {
            vect_[i] = &storage_[0] + i * ld_;
        }
    }

    void trade_boundaries()
    {
        assert(data_wghosts_ != 0);

        trade_data_tm_.start();

        data_wghosts_->trade_boundaries();

        trade_data_tm_.stop();
    }

    void scal(const int i, const double alpha);
    void scal(const double alpha);

    void set_zero();
    void set_zero(const int i, const short iloc);
    void set_zero(const int i);

    double dot(const int i, const int j, const short iloc) const;
    void scal(const double alpha, const int i, const short iloc);
    void axpy(const double alpha, const int ix, const int iy, const short iloc);
    void axpy(const double alpha, BlockVector& bv, const int ix, const int iy,
        const short iloc);

    BlockVector& operator-=(const BlockVector& src);

    void hasnan(const int j) const;

    int getld() const { return ld_; }

    static void setOverAllocateFactor(const float overallocate_factor)
    {
        assert(overallocate_factor > 0.);

        overallocate_factor_ = overallocate_factor;
    }

    static void incMaxAllocInstances(const short inc)
    {
        max_alloc_instances_ += inc;
        if (onpe0)
            std::cout << "BlockVector, max_alloc_instances_="
                      << max_alloc_instances_ << std::endl;

        // allocate extra memory now if initial blocks already allocated
        if (!class_storage_.empty())
            for (short i = 0; i < inc; ++i)
                allocate1NewBlock();
    }

    void set_ld_and_size_storage();

    static void printTimers(std::ostream& os);

    static int get_allocated_size_storage() { return allocated_size_storage_; }
};

template <typename ScalarType, typename MemorySpaceType>
std::vector<ScalarType*>
    BlockVector<ScalarType, MemorySpaceType>::class_storage_;
template <typename ScalarType, typename MemorySpaceType>
std::vector<short> BlockVector<ScalarType, MemorySpaceType>::allocated_;

// 4 slots for phi, residual, work, H*phi
template <typename ScalarType, typename MemorySpaceType>
short BlockVector<ScalarType, MemorySpaceType>::max_alloc_instances_ = 4;
template <typename ScalarType, typename MemorySpaceType>
pb::GridFuncVector<ScalarType, MemorySpaceType>*
    BlockVector<ScalarType, MemorySpaceType>::data_wghosts_
    = nullptr;

template <typename ScalarType, typename MemorySpaceType>
int BlockVector<ScalarType, MemorySpaceType>::size_storage_ = 0;
template <typename ScalarType, typename MemorySpaceType>
int BlockVector<ScalarType, MemorySpaceType>::numel_ = -1;
template <typename ScalarType, typename MemorySpaceType>
int BlockVector<ScalarType, MemorySpaceType>::locnumel_ = -1;
template <typename ScalarType, typename MemorySpaceType>
short BlockVector<ScalarType, MemorySpaceType>::subdivx_ = -1;
template <typename ScalarType, typename MemorySpaceType>
short BlockVector<ScalarType, MemorySpaceType>::n_instances_ = 0;
template <typename ScalarType, typename MemorySpaceType>
int BlockVector<ScalarType, MemorySpaceType>::allocated_size_storage_ = 0;
template <typename ScalarType, typename MemorySpaceType>
float BlockVector<ScalarType, MemorySpaceType>::overallocate_factor_
    = 0.; // should be explicitly set to value >= 1

template <typename ScalarType, typename MemorySpaceType>
int BlockVector<ScalarType, MemorySpaceType>::ld_ = 0;

template <typename ScalarType, typename MemorySpaceType>
Timer BlockVector<ScalarType, MemorySpaceType>::set_data_tm_(
    "BlockVector::set_data_wghosts");

template <typename ScalarType, typename MemorySpaceType>
Timer BlockVector<ScalarType, MemorySpaceType>::trade_data_tm_(
    "BlockVector::trade_data");

template <typename ScalarType, typename MemorySpaceType>
Timer BlockVector<ScalarType, MemorySpaceType>::assign_tm_(
    "BlockVector::assign");

template <typename ScalarType, typename MemorySpaceType>
Timer BlockVector<ScalarType, MemorySpaceType>::scal_tm_("BlockVector::scal");

template <typename ScalarType, typename MemorySpaceType>
Timer BlockVector<ScalarType, MemorySpaceType>::opminus_tm_(
    "BlockVector::opminus");

template <typename ScalarType, typename MemorySpaceType>
Timer BlockVector<ScalarType, MemorySpaceType>::copy_tm_("BlockVector::copy");

template <typename ScalarType, typename MemorySpaceType>
Timer BlockVector<ScalarType, MemorySpaceType>::diagop_tm_(
    "BlockVector::diagop");
#endif
