// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef BLOCKVECTOR_H
#define BLOCKVECTOR_H

#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "mputils.h"

#include "GridFuncVector.h"

#include <vector>

#define PAD     3

// memory allocation is done for the class and not for specific instances
// thus the use of static variables and functions
// this is done to avoid not finding blocks large enough later during run
// due to memory fragmentation
template <typename T>
class BlockVector
{
    static Timer  set_data_tm_;
    static Timer  trade_data_tm_;

    static short n_instances_;
    static short subdivx_;
    
    static pb::GridFuncVector<T>* data_wghosts_;
    
    // data allocator
    static std::vector<T*> class_storage_;
    
    // tells which blocks are allocated (1) or not (0) 
    // to a specific instantiated object
    static std::vector<short> allocated_;
    
    static short max_alloc_instances_;
    static int size_storage_;
    static int allocated_size_storage_;
    int size_storage_instance_;

    //safety coefficient to multiply dize of initial allocation
    //and have enough memory for later allocations
    static float overallocate_factor_;

    // index telling which slot of memory has been allocated 
    // for that object
    short my_allocation_;

    T* storage_;        // storage for data

    std::vector< T* > vect_;  // pointers to storage of vectors   
    
    static int numel_; // nb elements in each vector
    static int locnumel_;
    static int ld_;    // leading dimension
    
    int ld_instance_;
 
    const pb::Grid& mygrid_;
    short bc_[3];
    
    void allocate_storage();
    void deallocate_storage();
    BlockVector(const BlockVector& bv);
    void setup(const BlockVector& bv);
    
    static void allocate1NewBlock();
    
public:

    BlockVector(const pb::Grid& my_grid,
                const short subdivx,
                const short bc[3]);
    
    BlockVector(const BlockVector& bv, const bool copy_data);
    
    BlockVector& operator=(const BlockVector& bv);
    
    ~BlockVector();
    
    const pb::GridFuncVector<T>& getDataWGhosts() 
    {
        assert( data_wghosts_!=0 );
        return *data_wghosts_;
    }
    
    pb::GridFuncVector<T>* getPtDataWGhosts() 
    {
        return data_wghosts_;
    }
    
    void initialize(const std::vector<std::vector<int> >& gid,
                    const bool skinny_stencil);
    
    void clear()
    {
        vect_.clear();
        
        deallocate_storage();
    }
    
    void axpy(const double alpha, const BlockVector& bv)
    {
        assert( storage_!=0 );
        assert( bv.storage_!=0 );
        assert( storage_!=bv.storage_ );
    
        MPaxpy(size_storage_, alpha, bv.storage_, storage_);
    }
    
    void axpy(const double alpha, const T* const vy)
    {
        MPaxpy(size_storage_, alpha, vy, storage_);
    }
    
    T* vect(const int i)const
    {
        assert( i<(int)vect_.size() );
        assert( vect_[i]!=0 );
        return vect_[i];
    }
    
    T maxAbsValue()const;

    template<typename T2>
    void setDataWithGhosts(pb::GridFuncVector<T2>* data_wghosts);

    void setDataWithGhosts();
    
    void assign(const int color, const T* const src, const int n=1)
    {
        assert( (color+n-1) < (int)vect_.size() );
        int ione=1;
        int mysize=n*ld_;
        Tcopy(&mysize, src, &ione, vect_[color], &ione);
    }
    
    void assign(const T* const src)
    {
        memcpy(storage_, src, size_storage_*sizeof(T));
    }

    /*
     * assign functions for source data with ghost values
     */
    template<typename T2>
    void assign(const pb::GridFuncVector<T2>& src);
    template<typename T2>
    void assignComponent(const pb::GridFunc<T2>& src, const int i);
    template<typename T2>
    void assignComponent(const pb::GridFuncVector<T2>& src, const int i);

    void assignLocal(const int color, const short iloc,
                     const T* const src)
    {
        assert( color>=0 );
        assert( color<(int)vect_.size() );
        assert( iloc<subdivx_ );
        memcpy(vect_[color]+iloc*locnumel_, 
               src, locnumel_*sizeof(T));
    }
    
    void copyDataFrom(const BlockVector& src)
    {
        assert( src.size_storage_==size_storage_ );
        assert( storage_!=0 );
        assert( src.storage_!=0 );
        memcpy(storage_, src.storage_, size_storage_*sizeof(T));
    }
    
    pb::GridFunc<T>& getVectorWithGhosts(const int i)
    {
        assert( data_wghosts_!=0 );
        assert( i<vect_.size() );
    
        return data_wghosts_->func(i);
    }

    void setStorage(T* new_storage)
    {
        assert( new_storage!=0 || vect_.size()==0 );
        
        storage_=new_storage;
        for(int i = 0;i < vect_.size();i++)
        {
            vect_[i] = &storage_[0] + i*ld_;
        }   
    }
    
    void trade_boundaries()
    {
        assert( data_wghosts_!=0 );
        
        trade_data_tm_.start();
        
        //if( onpe0 )
        //    (*MPIdata::sout)<<"BlockVector::trade_boundaries()"<<endl;
        data_wghosts_->trade_boundaries();
        
        trade_data_tm_.stop();
    }

    void scal(const int i, const double alpha);
    void scal(const double alpha);
    void set_zero()
    {
        memset(storage_, 0, size_storage_*sizeof(T));
    }
    void set_zero(const int i, const short iloc)
    {
        assert( i<vect_.size() );
        assert( iloc<subdivx_ );
        memset(vect_[i]+iloc*locnumel_, 0, locnumel_*sizeof(T));
    }
    void set_zero(const int i)
    {
        assert( i<vect_.size() );
        
        for(short iloc=0;iloc<subdivx_;iloc++)
        memset(vect_[i]+iloc*locnumel_, 0, locnumel_*sizeof(T));
    }
    double dot(const int i, const int j, const short iloc)const;
    void scal(const double alpha, const int i, const short iloc);
    void axpy(const double alpha,const int ix, const int iy, const short iloc);
    void axpy(const double alpha, BlockVector& bv, const int ix, const int iy, const short iloc);
    
    BlockVector& operator-=(const BlockVector& src);

    void hasnan(const int j)const;
    
    int getld()const{ return ld_; }
    
    static void setOverAllocateFactor(const float overallocate_factor)
    {
        assert( overallocate_factor>0. );
        
        overallocate_factor_=overallocate_factor;
    }
    
    static void incMaxAllocInstances(const short inc)
    {
        max_alloc_instances_+=inc;
        if( onpe0 )std::cout<<"BlockVector, max_alloc_instances_="
                            <<max_alloc_instances_<<std::endl;
        
        // allocate extra memory now if initial blocks already allocated
        if( !class_storage_.empty() )
            for(short i=0;i<inc;++i)allocate1NewBlock();
    }

    void set_ld_and_size_storage();

    static void printTimers(std::ostream& os);
};

template <typename T>
std::vector<T*>       BlockVector<T>::class_storage_;
template <typename T>
std::vector<short> BlockVector<T>::allocated_;
template <typename T>
short         BlockVector<T>::max_alloc_instances_=7;
template <typename T>
pb::GridFuncVector<T>* BlockVector<T>::data_wghosts_=0;

template <typename T>
int BlockVector<T>::size_storage_=0;
template <typename T>
int BlockVector<T>::numel_=-1;
template <typename T>
int BlockVector<T>::locnumel_=-1;
template <typename T>
short BlockVector<T>::subdivx_=-1;
template <typename T>
short BlockVector<T>::n_instances_=0;
template <typename T>
int BlockVector<T>::allocated_size_storage_=0;
template <typename T>
float BlockVector<T>::overallocate_factor_=0.; // should be explicitly set to value >= 1

template <typename T>
int BlockVector<T>::ld_=0;

template <typename T>
Timer BlockVector<T>::set_data_tm_("BlockVector::set_data_wghosts");

template <typename T>
Timer BlockVector<T>::trade_data_tm_("BlockVector::trade_data");
#endif
