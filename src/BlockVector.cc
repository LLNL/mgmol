// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "BlockVector.h"
#include "Control.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "global.h"
#include "memory_space.h"

#define NEWSTORAGE 0

template <typename ScalarType, typename MemorySpaceType>
BlockVector<ScalarType, MemorySpaceType>::BlockVector(
    const pb::Grid& my_grid, const short subdivx, const short bc[3])
    : mygrid_(my_grid)
{
    for (short i = 0; i < 3; i++)
        bc_[i] = bc[i];
    if (n_instances_ == 0)
    {
        subdivx_  = subdivx;
        numel_    = (int)my_grid.size();
        locnumel_ = numel_ / subdivx_;
        ld_       = numel_ + PAD;
    }

    ld_instance_ = ld_;

    storage_       = nullptr;
    my_allocation_ = -1;

    n_instances_++;
}
template <typename ScalarType, typename MemorySpaceType>
BlockVector<ScalarType, MemorySpaceType>::~BlockVector()
{
    deallocate_storage();
    if (n_instances_ == 1)
    {
        delete data_wghosts_;
        data_wghosts_ = nullptr;
        allocated_.clear();
        for (typename std::vector<ScalarType*>::iterator it
             = class_storage_.begin();
             it != class_storage_.end(); ++it)
            delete[] * it;
    }
    n_instances_--;
}

template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::allocate1NewBlock()
{
    assert(size_storage_ > 0);
    assert(overallocate_factor_ >= 1.);

    // size of allocation is based on first call to this function
    if (allocated_size_storage_ == 0)
        allocated_size_storage_
            = (int)((ScalarType)(size_storage_)*overallocate_factor_);

    // check if storage required is not bigger than initial allocation
    if (size_storage_ > allocated_size_storage_)
    {
        (*MPIdata::sout)
            << "BlockVector::allocate1NewBlock(): allocated_size_storage_ = "
            << allocated_size_storage_ << std::endl;
        (*MPIdata::sout)
            << "BlockVector::allocate1NewBlock(): size_storage_           = "
            << size_storage_ << std::endl;
        (*MPIdata::sout)
            << "BlockVector::allocate1NewBlock(): size_storage_ too large!!!"
            << std::endl;
        exit(0);
    }

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "BlockVector, allocate memory size "
                         << allocated_size_storage_ << std::endl;

    class_storage_.push_back(new ScalarType[allocated_size_storage_]);

    memset((class_storage_.back()), 0,
        allocated_size_storage_ * sizeof(ScalarType));
}
template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::allocate_storage()
{
    assert(size_storage_ > 0);
    assert(storage_ == nullptr);

    Control& ct                  = *(Control::instance());
    static short high_water_mark = 0;

#ifdef NEWSTORAGE
    // firt time, allocate memory for max_alloc_instances_ objects
    if (allocated_.size() == 0)
    {
        if (onpe0 && ct.verbose > 1)
            (*MPIdata::sout) << "BlockVector, allocate " << max_alloc_instances_
                             << " instances" << std::endl;
        for (short i = 0; i < max_alloc_instances_; i++)
            allocate1NewBlock();

        allocated_.clear();
        allocated_.push_back(1);
        my_allocation_ = 0;
    }
    else
    {
        std::vector<short>::iterator ia = allocated_.begin();
        short ii                        = 0;
        while (ia != allocated_.end())
        {
            if (*ia == 0) break;
            ia++;
            ii++;
        }
        if (ia == allocated_.end())
        {
            allocated_.push_back(1);
        }
        else
        {
            *ia = 1;
        }
        my_allocation_ = ii;
        if (my_allocation_ >= static_cast<int>(class_storage_.size()))
        {
            std::cerr << "Number of allocated blocks: " << class_storage_.size()
                      << ", Cannot allocate another BlockVector!!!"
                      << std::endl;
            exit(0);
        }
    }

    if (my_allocation_ > high_water_mark) high_water_mark = my_allocation_;
    if (onpe0 && ct.verbose > 2)
    {
        (*MPIdata::sout) << "allocation in slot " << my_allocation_
                         << std::endl;
        (*MPIdata::sout) << "high water mark: " << high_water_mark << std::endl;
    }
    assert(my_allocation_ < static_cast<int>(class_storage_.size()));
    assert(my_allocation_ < static_cast<int>(class_storage_.size()));
    assert(my_allocation_ >= 0);

    if (allocated_size_storage_ < size_storage_)
    {
        std::cerr << "ERROR BlockVector: trying to use allocation "
                  << size_storage_ << " bigger than initialy preallocated "
                  << allocated_size_storage_ << "!!!" << std::endl;
        ct.global_exit(0);
    }
    storage_ = class_storage_[my_allocation_];
    assert(storage_ != nullptr);
#else
    storage_ = new ScalarType[size_storage_];
#endif

    assert(storage_ != nullptr);
    assert(class_storage_.size() > 0);
}

// free slot of memory allocated to particular object
template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::deallocate_storage()
{
    if (my_allocation_ >= 0)
    {
        assert(storage_ != nullptr);
#ifdef NEWSTORAGE
        assert(my_allocation_ >= 0);
        allocated_[my_allocation_] = 0;
        my_allocation_             = -1;
        storage_                   = nullptr;
#else
        delete[] storage_;
#endif
        storage_ = nullptr;
    }
}
template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::setup(const BlockVector& bv)
{
    storage_       = nullptr;
    my_allocation_ = -1;
    for (short i = 0; i < 3; i++)
        bc_[i] = bv.bc_[i];
    allocate_storage();

    // set vect_
    vect_.resize(bv.vect_.size());
    for (unsigned int i = 0; i < vect_.size(); i++)
    {
        assert(i * ld_ < static_cast<unsigned int>(size_storage_));
        vect_[i] = storage_ + i * ld_;
    }
}
template <typename ScalarType, typename MemorySpaceType>
BlockVector<ScalarType, MemorySpaceType>::BlockVector(
    const BlockVector<ScalarType, MemorySpaceType>& bv, const bool copy_data)
    : size_storage_instance_(bv.size_storage_instance_),
      ld_instance_(bv.ld_instance_),
      mygrid_(bv.mygrid_)
{
    n_instances_++;

    set_ld_and_size_storage();

    setup(bv);

    if (copy_data)
        memcpy(storage_, bv.storage_, size_storage_ * sizeof(ScalarType));
}
template <typename ScalarType, typename MemorySpaceType>
BlockVector<ScalarType, MemorySpaceType>&
BlockVector<ScalarType, MemorySpaceType>::operator=(
    const BlockVector<ScalarType, MemorySpaceType>& bv)
{
    if (this == &bv) return *this;

    n_instances_++;

    setup(bv);

    memcpy(storage_, bv.storage_, size_storage_ * sizeof(ScalarType));

    return *this;
}
template <typename ScalarType, typename MemorySpaceType>
BlockVector<ScalarType, MemorySpaceType>&
BlockVector<ScalarType, MemorySpaceType>::operator-=(
    const BlockVector<ScalarType, MemorySpaceType>& src)
{
    for (unsigned int i = 0; i < vect_.size(); i++)
    {
        ScalarType* vi             = vect_[i];
        const ScalarType* const si = src.vect_[i];
        for (int j = 0; j < numel_; j++)
            vi[j] -= si[j];
    }
    return *this;
}
template <typename ScalarType, typename MemorySpaceType>
template <typename ScalarType2>
void BlockVector<ScalarType, MemorySpaceType>::assign(
    const pb::GridFuncVector<ScalarType2>& src)
{
    for (unsigned int i = 0; i < vect_.size(); i++)
    {
        ScalarType* dest = vect_[i];
        src.getValues(i, dest);
    }
}

template <typename ScalarType, typename MemorySpaceType>
template <typename ScalarType2>
void BlockVector<ScalarType, MemorySpaceType>::assignComponent(
    const pb::GridFunc<ScalarType2>& src, const int i)
{
    ScalarType* dest = vect_[i];
    src.getValues(dest);
}

template <typename ScalarType, typename MemorySpaceType>
template <typename ScalarType2>
void BlockVector<ScalarType, MemorySpaceType>::assignComponent(
    const pb::GridFuncVector<ScalarType2>& src, const int i)
{
    ScalarType* dest = vect_[i];
    src.getValues(i, dest);
}

template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::initialize(
    const std::vector<std::vector<int>>& gid, const bool skinny_stencil)
{
    assert(storage_ == nullptr);

    int nbvect = (int)gid[0].size();

    // set minimum number of vectors to take care of case of empty subdomains
    const short nvec       = nbvect > 20 ? nbvect : 20;
    size_storage_          = ld_ * nvec;
    size_storage_instance_ = size_storage_;

    allocate_storage();
    memset(storage_, 0, size_storage_ * sizeof(ScalarType));

    vect_.resize(nbvect);
    for (int i = 0; i < nbvect; i++)
    {
        vect_[i] = storage_ + i * ld_;
    }

    // rebuild data_wghosts_ with new gid
    if (data_wghosts_ != nullptr)
    {
        delete data_wghosts_;
        data_wghosts_ = nullptr;
    }

    data_wghosts_ = new pb::GridFuncVector<ScalarType>(
        true, mygrid_, bc_[0], bc_[1], bc_[2], gid, skinny_stencil);

    data_wghosts_->resetData();

    data_wghosts_->set_updated_boundaries(false);
}
template <typename ScalarType, typename MemorySpaceType>
double BlockVector<ScalarType, MemorySpaceType>::dot(
    const int i, const int j, const short iloc) const
{
    assert(i < static_cast<int>(vect_.size()));
    assert(j < static_cast<int>(vect_.size()));
    assert(iloc < subdivx_);
    assert(vect_[i] != nullptr);
    assert(vect_[j] != nullptr);

    const int shift = iloc * locnumel_;
    return MPdot(locnumel_, vect_[i] + shift, vect_[j] + shift);
}
template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::scal(
    const int i, const double alpha)
{
    assert(i < static_cast<int>(vect_.size()));
    assert(vect_[i] != nullptr);

    MPscal(ld_, alpha, vect_[i]);
}
template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::scal(const double alpha)
{
    assert(storage_ != nullptr);

    MPscal(size_storage_, alpha, storage_);
}
template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::scal(
    const double alpha, const int i, const short iloc)
{
    assert(static_cast<unsigned int>(i) < vect_.size());
    assert(vect_[i] != nullptr);

    const int shift = iloc * locnumel_;
    MPscal(locnumel_, alpha, vect_[i] + shift);
}
template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::axpy(
    const double alpha, const int ix, const int iy, const short iloc)
{
    assert(static_cast<unsigned int>(ix) < vect_.size());
    assert(static_cast<unsigned int>(iy) < vect_.size());
    assert(vect_[ix] != nullptr);
    assert(vect_[iy] != nullptr);
    assert(vect_[ix] != vect_[iy]);

    const int shift = iloc * locnumel_;

    MPaxpy(locnumel_, alpha, vect_[ix] + shift, vect_[iy] + shift);
}
template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::axpy(const double alpha,
    BlockVector<ScalarType, MemorySpaceType>& bv, const int ix, const int iy,
    const short iloc)
{
    assert(ix < (int)bv.vect_.size());
    assert(iy < (int)vect_.size());
    assert(bv.vect_[ix] != vect_[iy]);

    const int shift = iloc * locnumel_;

    MPaxpy(locnumel_, alpha, bv.vect_[ix] + shift, vect_[iy] + shift);
}
template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::hasnan(const int j) const
{
    for (int i = 0; i < ld_; i++)
    {
        if (std::isnan(vect_[j][i]))
            (*MPIdata::sout) << "BlockVector: Nan in column " << j << ", row "
                             << i << std::endl;
    }
}
template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::setDataWithGhosts()
{
    // set internal data
    setDataWithGhosts(data_wghosts_);
}

template <typename ScalarType, typename MemorySpaceType>
template <typename ScalarType2>
void BlockVector<ScalarType, MemorySpaceType>::setDataWithGhosts(
    pb::GridFuncVector<ScalarType2>* data_wghosts)
{
    assert(data_wghosts != nullptr);

    set_data_tm_.start();

    data_wghosts->resetData();

    data_wghosts->set_updated_boundaries(false);

    for (int i = 0; i < (int)vect_.size(); i++)
    {
        data_wghosts->assign(i, vect_[i]);
    }

    set_data_tm_.stop();
}

template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::printTimers(std::ostream& os)
{
    set_data_tm_.print(os);
    trade_data_tm_.print(os);
}

template <typename ScalarType, typename MemorySpaceType>
ScalarType BlockVector<ScalarType, MemorySpaceType>::maxAbsValue() const
{
    int ione        = 1;
    int imax        = IDAMAX(&size_storage_, storage_, &ione);
    ScalarType maxv = fabs(storage_[imax - 1]);
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&maxv, 1, MPI_MAX);
    return maxv;
}

template <typename ScalarType, typename MemorySpaceType>
void BlockVector<ScalarType, MemorySpaceType>::set_ld_and_size_storage()
{
    size_storage_ = size_storage_instance_;
    ld_           = ld_instance_;
}

template class BlockVector<double, MemorySpace::Host>;
template void BlockVector<double, MemorySpace::Host>::assign(
    const pb::GridFuncVector<float>& src);
template void BlockVector<double, MemorySpace::Host>::assign(
    const pb::GridFuncVector<double>& src);
template void BlockVector<double, MemorySpace::Host>::assignComponent(
    const pb::GridFunc<float>& src, const int i);
template void BlockVector<double, MemorySpace::Host>::assignComponent(
    const pb::GridFunc<double>& src, const int i);
template void BlockVector<double, MemorySpace::Host>::setDataWithGhosts(
    pb::GridFuncVector<float>* data_wghosts);
template void BlockVector<double, MemorySpace::Host>::setDataWithGhosts(
    pb::GridFuncVector<double>* data_wghosts);
#ifdef USE_MP
template class BlockVector<float, MemorySpace::Host>;
template void BlockVector<float, MemorySpace::Host>::assign(
    const pb::GridFuncVector<float>& src);
template void BlockVector<float, MemorySpace::Host>::assign(
    const pb::GridFuncVector<double>& src);
template void BlockVector<float, MemorySpace::Host>::assignComponent(
    const pb::GridFunc<float>& src, const int i);
template void BlockVector<float, MemorySpace::Host>::assignComponent(
    const pb::GridFunc<double>& src, const int i);
template void BlockVector<float, MemorySpace::Host>::setDataWithGhosts(
    pb::GridFuncVector<float>* data_wghosts);
template void BlockVector<float, MemorySpace::Host>::setDataWithGhosts(
    pb::GridFuncVector<double>* data_wghosts);
#endif
