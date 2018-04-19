// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef KBPROJECTOR_H
#define KBPROJECTOR_H

#include "Mesh.h"
#include "mputils.h"
#include "GridFunc.h"
#include "LocGridOrbitals.h"
#include "global.h"
#include <vector>
#include <cassert>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

class Species;

#define CHECK_NORM  0

// one KB projector
class KBprojectorSparse
{
    // work arrays (1 for each thread)
    static std::vector< std::vector<ORBDTYPE> > work_nlindex_;

    static std::vector< std::vector<KBPROJDTYPE> > work_proj_;

    // reference to species associated to projector
    const Species& species_;
    
    short subdivx_;
    short maxl_;
    short llocal_;
    
    // projectors for iloc, l, p, m
    std::vector< std::vector< std::vector< std::vector< std::vector<KBPROJDTYPE> > > > > projector_;
    
    std::vector<int> size_nl_;
    
    short range_kbproj_;

    double center_[3];
    
    // projector start
    short kb_proj_start_index_[3];
    double kb_proj_start_[3];
    
    std::vector<short> proj_indices_[3];
    
    // list of nodes where the non-local projector is non-zero on this PE
    // for each subdomain "iloc"
    std::vector< std::vector<int> >    nlindex_;
    
    bool** is_in_domain_;
    
    // multiplicity of projector for each "l"
    std::vector<short> multiplicity_;
    
#if CHECK_NORM
    std::vector< std::vector< double > > norm2_;
#endif
    
    // private functions
    void setProjectors(const short iloc, const int icount);
    
    void setSProjector(const short iloc, const int icount);
    void setPProjector(const short iloc, const int icount);
    void setDProjector(const short iloc, const int icount);
    void setFProjector(const short iloc, const int icount);
    
    void setNLindex(const short iloc, const int size,
                    const int* const pvec);

    void setKBProjStart();
    void setProjIndices(const short dir);
    int get_index_array(int* pvec, const short iloc,
                        const short index_low[3], 
                        const short index_high[3]);
    bool overlapWithBox(const short index_low[3], 
                        const short index_high[3]);

    //get projector l,m in pieces corresponding to subdivisions
    const KBPROJDTYPE* getProjector(const short iloc, const short l, const short p, const short m)const
    {
        assert( iloc<subdivx_ );
        assert( l<projector_[iloc].size() );
        assert( p<projector_[iloc][l].size() );
        assert( m<projector_[iloc][l][p].size() );
        
        return &projector_[iloc][l][p][m][0];
    }
    std::vector<KBPROJDTYPE>& getRefProjector(const short iloc, const short l, const short p, const short m)
    {
        assert( iloc<subdivx_ );
        assert( l<projector_[iloc].size() );
        assert( p<projector_[iloc][l].size() );
        assert( m<projector_[iloc][l][p].size() );
        
        return projector_[iloc][l][p][m];
    }
    double dotPsi(const short iloc, const short l, const short p, const short m)const
    {
        assert( iloc<subdivx_ );
        assert( l<projector_[iloc].size() );
        assert( p<projector_[iloc][l].size() );
        assert( m<projector_[iloc][l][p].size() );
        assert( omp_get_thread_num()<work_nlindex_.size() );
        
        return MPdot(size_nl_[iloc], &work_nlindex_[omp_get_thread_num()][0], 
                                     &projector_[iloc][l][p][m][0]);
    }
    bool setIndexesAndProjectors();
    void allocateProjectors();

    void getProjectors(const short iloc, 
                       std::vector<const KBPROJDTYPE*>& projectors)const;

public:
    KBprojectorSparse(const Species& sp);
    KBprojectorSparse(const KBprojectorSparse& kbp);

    ~KBprojectorSparse()
    {
        clear();
    }
    
    void clear()
    {
        if( is_in_domain_!=NULL )
        {
            for(short iloc=0;iloc<subdivx_;iloc++)
                delete[] is_in_domain_[iloc];
            delete[] is_in_domain_; is_in_domain_=NULL;
        }
        
        for(short dir=0;dir<3;dir++)proj_indices_[dir].clear();
    }
    void setup(const short subdivx, const double center[3]);
    void initCenter(const double center[3])
    {
        for(short i=0;i<3;i++)
            center_[i]=center[i];
    }
    
    double maxRadius()const
    {
        Mesh* mymesh = Mesh::instance();
        const pb::Grid& mygrid  = mymesh->grid();

        double radius=0.;
        
        if( maxl_>0 )
        for(short dir=0;dir<3;dir++)
        {
            const double hgrid = mygrid.hgrid(dir);
        
            double left =center_[dir]-kb_proj_start_[dir];
            double right=kb_proj_start_[dir]+(range_kbproj_-1)*hgrid-center_[dir];
            
            radius = left  > radius ? left  : radius;
            radius = right > radius ? right : radius;
            
            //cout<<"kb_proj_start_[dir]="<<kb_proj_start_[dir]
            //    <<", center="<<center_[dir]
            //    <<", left="<<left<<", right="<<right<<endl;
        }
        
        return radius;
    }
    
    bool overlapPE()const;
    void init_work_nlindex(const short iloc, const ORBDTYPE* const psi);
    
    bool overlaps(const short iloc)const
    {
        return (size_nl_[iloc]>0);
    }
    double dotPsi(const short iloc, const short index)const;
    
    //axpySket for templated destination type
    template <typename T>
    void axpySKet(const short iloc, 
                  const double alpha, T* const)const;
    template <typename T>
    void axpyKet(const short iloc, 
                 const vector<double>& alpha, 
                 T* const dst)const;
    
    bool onlyOneProjector()const
    {
        return ( (maxl_ == 1) && (llocal_==1) && (multiplicity_[0]==1) );
    }
   
    short nProjectors()const
    {
        short nproj=0;
        assert(llocal_<4);
        for(short l=0;l<=maxl_;l++)
            if(llocal_!=l)
                nproj+=(2*l+1)*multiplicity_[l];
        return nproj;
    }
       
    short nProjectorsSubdomain()const
    {
        short nproj=0;
        assert(llocal_<4);
        if( overlapPE() )
            nproj=nProjectors();
        return nproj;
    }
       
    void getKBsigns(std::vector<short>& kbsigns)const;
    void getKBcoeffs(std::vector<double>& coeffs)const;
};

#endif
