// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: GridFunc<T>.cc,v 1.42 2010/01/28 22:56:31 jeanluc Exp $
#include <cstdlib>
#include <cstdarg>
#include <fstream>
#include "GridFunc.h"
#include "radial_functions.h"
#include "MGmol_blas1.h"
#include "mputils.h"
#include "MGmol_MPI.h"

#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

const double inv64=1./64.;

#ifdef EXPLICIT_INSTANTIATION_VINT
namespace std{
template class vector<int>;
template class vector<vector<int> >;
}
#endif

namespace pb{

Timer  GridFuncInterface::trade_bc_tm_("GridFunc::trade_bc");
Timer  GridFuncInterface::restrict3D_tm_("GridFunc::restrict3D");
Timer  GridFuncInterface::extend3D_tm_("GridFunc::extend3D");
Timer  GridFuncInterface::prod_tm_("GridFunc::prod");
Timer  GridFuncInterface::gather_tm_("GridFunc::gather");
Timer  GridFuncInterface::scatter_tm_("GridFunc::scatter");
Timer  GridFuncInterface::all_gather_tm_("GridFunc::all_gather");
Timer  GridFuncInterface::finishExchangeNorthSouth_tm_("GridFunc::finishExNorthSouth");
Timer  GridFuncInterface::finishExchangeUpDown_tm_("GridFunc::finishExUpDown");
Timer  GridFuncInterface::finishExchangeEastWest_tm_("GridFunc::finishExEastWest");

template <typename T>vector<T> GridFunc<T>::buf1_;
template <typename T>vector<T> GridFunc<T>::buf2_;
template <typename T>vector<T> GridFunc<T>::buf3_;
template <typename T>vector<T> GridFunc<T>::buf4_;

template <typename T>
T GridFunc<T>::ran0()
{
    static bool first_time=true;
    if( first_time )
    {
        srand(1);
        first_time=false;
    }
    
    return static_cast<T>(std::rand())/static_cast<T>(RAND_MAX);
}

template <typename T>
T GridFunc<T>::radial_func(const double r, const double a, const short ftype)
{
    if(ftype == 0) return one(r, a);
    //add other options here
    else return 1.0;
}

template <>
double GridFunc<double>::fmax()
{
    int ione=1;
    int n=grid_.sizeg();
    int imax=idamax(&n, &uu_[0], &ione)-1;
    
    double vmax=mype_env().double_max_all(uu_[imax]);

    return vmax;
}

template <>
double GridFunc<float>::fmax()
{
    int ione=1;
    int n=grid_.sizeg();
    int imax=isamax(&n, &uu_[0], &ione)-1;
    
    double vmax=mype_env().double_max_all((double)uu_[imax]);

    return vmax;
}

template <typename T>
void GridFunc<T>::setup()
{
    directionDirichlet_[0]=(bc_[0]==0);
    directionDirichlet_[1]=(bc_[1]==0);
    directionDirichlet_[2]=(bc_[2]==0);

    directionPeriodic_[0] =(bc_[0]==1);
    directionPeriodic_[1] =(bc_[1]==1);
    directionPeriodic_[2] =(bc_[2]==1);

    directionMultipole_[0]=(bc_[0]==2);
    directionMultipole_[1]=(bc_[1]==2);
    directionMultipole_[2]=(bc_[2]==2);

    directionNeumann_[0]  =(bc_[0]==3);
    directionNeumann_[1]  =(bc_[1]==3);
    directionNeumann_[2]  =(bc_[2]==3);
    
    mytask_=mype_env().mytask();

    dim_[0] = grid_.dim(0);
    dim_[1] = grid_.dim(1);
    dim_[2] = grid_.dim(2);
    
    incx_=grid_.inc(0);
    incy_=grid_.inc(1);

    //set flags to true if
    //periodic in y direction or
    //or if there is another subdomain to get data from in that direction
    north_=(directionPeriodic_[1] || (mype_env().mpi_neighbors(NORTH)>mytask_));
    south_=(directionPeriodic_[1] || (mype_env().mpi_neighbors(SOUTH)<mytask_));

    up_  =( directionPeriodic_[2] || (mype_env().mpi_neighbors(UP)>mytask_) );
    down_=( directionPeriodic_[2] || (mype_env().mpi_neighbors(DOWN)<mytask_) );

    east_ =( directionPeriodic_[0] || (mype_env().mpi_neighbors(EAST)>mytask_) );
    west_ =( directionPeriodic_[0] || (mype_env().mpi_neighbors(WEST)<mytask_) );

    assert(dim_[0]>=ghost_pt());
    assert(dim_[1]>=ghost_pt());
    assert(dim_[2]>=ghost_pt());    
    assert(grid_.inc(2)==1);
}

template <typename T>
GridFunc<T>::GridFunc(const Grid& my_grid,
                      const short px, const short py, const short pz):
          grid_(my_grid)
{
    assert( px==0 || px==1 || px==2 ); 
    assert( py==0 || py==1 || py==2 ); 
    assert( pz==0 || pz==1 || pz==2 ); 

    bc_[0]=px;
    bc_[1]=py;
    bc_[2]=pz;

    setup();

    for(short i=0;i<3;i++){
        assert( grid_.dim(i)>0 );
        assert( grid_.dim(i)<10000 );
    }        
    
    //cout<<"Constructor with BC for GridFunc<T> at level "<<grid_.level()<<endl;
    updated_boundaries_=true; // boundaries initialized

    //cout<<"nn="<<grid_.sizeg<<endl;
    if( grid_.active() )
    {
        const size_t n=grid_.sizeg();
        uu_ = new T[n];
        memset(uu_,0,n*sizeof(T));
    }else{
        uu_=0;
    }
    const short shift=ghost_pt();
    const int dimxy=(dim(1)+2*ghost_pt())*dim(0);
    
    // resize static buffers if needed
    if(mype_env().n_mpi_tasks()>1){
        int size_max = shift*grid_.inc(0);
        size_max = max(size_max, shift*dimxy);
        size_max = max(size_max, shift*dim(0)*grid_.inc(1));

        if( size_max>(int)buf1_.size() ){
            buf1_.resize(size_max);
            buf2_.resize(size_max);
            buf3_.resize(size_max);
            buf4_.resize(size_max);
        }
    }
}

// copy constructor - double
template <typename T>
GridFunc<T>::GridFunc(const GridFunc<double>& A):grid_(A.grid())
{
    assert( grid_.dim(0)>0 );
    assert( grid_.dim(0)<10000 );

    bc_[0]=A.bc(0);
    bc_[1]=A.bc(1);
    bc_[2]=A.bc(2);

    setup();

    if( !grid_.active() )return;
    int n=grid_.sizeg();
    
    double * au = A.uu();
    if(au != 0){
        uu_ = new T[n];
        MPcpy(uu_, au, n);
        updated_boundaries_=A.updated_boundaries();
    }else{
        updated_boundaries_=false;
        uu_ = NULL;
    }

    //cout<<"Copy constructor for function on grid "<<grid_.level()<<endl;    

}

// copy constructor - float
template <typename T>
GridFunc<T>::GridFunc(const GridFunc<float>& A):grid_(A.grid())
{
    assert( grid_.dim(0)>0 );
    assert( grid_.dim(0)<10000 );

    bc_[0]=A.bc(0);
    bc_[1]=A.bc(1);
    bc_[2]=A.bc(2);

    setup();

    if( !grid_.active() )return;
    int n=grid_.sizeg();
    
    float * au = A.uu();
    if(au != 0){
        uu_ = new T[n];
        MPcpy(uu_, au, n);
        updated_boundaries_=A.updated_boundaries();
    }else{
        updated_boundaries_=false;
        uu_ = NULL;
    }

    //cout<<"Copy constructor for function on grid "<<grid_.level()<<endl;    

}

// copy constructor on different grid
template <typename T>
GridFunc<T>::GridFunc(const GridFunc<T>& A, const Grid& new_grid):grid_(new_grid)
{
    assert(dim(0)>=1);
    
    bc_[0]=A.bc_[0];
    bc_[1]=A.bc_[1];
    bc_[2]=A.bc_[2];
    
    setup();

    if( !grid_.active() )return;

    T  *vv=A.uu_;
    const short shift1=ghost_pt();
    const short shift2=A.ghost_pt();

    const int incx2=A.grid_.inc(0);
    const int incy2=A.grid_.inc(1);

    assert(grid_.inc(2)==1);
    assert(A.grid_.inc(2)==1);
    
    if(A.uu_!=0){
        uu_ = new T[grid_.sizeg()];
        memset(uu_,0,grid_.sizeg()*sizeof(T));
        size_t sdim2=dim_[2]*sizeof(T);
    
        for(int ix = 0;ix < dim_[0];ix++) {

            int ix1 = (ix+ shift1)*incx_+shift1;
            int ix2 = (ix+ shift2)*incx2+shift2;

            for(int iy = 0;iy < dim_[1];iy++) {

                int iy1 = ix1+(iy+ shift1)*incy_;
                int iy2 = ix2+(iy+ shift2)*incy2;
    
                memcpy(uu_+iy1, vv+iy2, sdim2);
            } 

        } 

    
        updated_boundaries_=false;
    }else{
        cout<<"Invalid initialization"<<endl;
        exit(2);
    }

    //cout<<"Copy constructor for function on grid "<<grid_.level()<<endl;    

}

// Constructor
template <typename T>
GridFunc<T>::GridFunc(const T* const vv, const Grid& new_grid, 
                   const short px, const short py, const short pz):grid_(new_grid)
{
    assert( px==0 || px==1 || px==2 ); 
    assert( py==0 || py==1 || py==2 ); 
    assert( pz==0 || pz==1 || pz==2 ); 

    bc_[0]=px;
    bc_[1]=py;
    bc_[2]=pz;
    
    setup();

    if( !grid_.active() )return;

    const short nghosts=ghost_pt();

    const int incx2=dim(2)*dim(1);
    const int incy2=dim(2);

    assert(grid_.inc(2)==1);
 
    const size_t sdim2=dim_[2]*sizeof(T);
    if(vv!=0){
        uu_ = new T[grid_.sizeg()];
        memset(uu_,0,grid_.sizeg()*sizeof(T));

        for(int ix = 0;ix < dim_[0];ix++) {

            const int ix1 = (ix+ nghosts)*incx_
                          + nghosts*(1+incy_);
            const int ix2 = ix*incx2;

            for(int iy = 0;iy < dim_[1];iy++) {

                int iy1 = ix1+iy*incy_;
                int iy2 = ix2+iy*incy2;
    
                memcpy(uu_+iy1, vv+iy2, sdim2);
            } 

        }

        updated_boundaries_=false;
    }else{
        updated_boundaries_=false;
        uu_ = 0;
        //cout<<" Empty GridFunc<T>!!!"<<endl;
    }

    //cout<<"Copy constructor for function on grid "<<grid_.level()<<endl;    

    const short shift=ghost_pt();
    const int dimxy=(dim(1)+2*ghost_pt())*dim(0);
    
    // resize static buffers if needed
    if(mype_env().n_mpi_tasks()>1){
        int size_max = shift*grid_.inc(0);
        size_max = max(size_max, shift*dimxy);
        size_max = max(size_max, shift*dim(0)*grid_.inc(1));

        if( size_max>(int)buf1_.size() ){
            buf1_.resize(size_max);
            buf2_.resize(size_max);
            buf3_.resize(size_max);
            buf4_.resize(size_max);
        }
    }
}

// Constructor
template <typename T>
GridFunc<T>::GridFunc(const T* const src, const Grid& new_grid, 
                   const short px, const short py, const short pz,
                   const char dis):grid_(new_grid)
{
    assert( px==0 || px==1 || px==2 ); 
    assert( py==0 || py==1 || py==2 ); 
    assert( pz==0 || pz==1 || pz==2 ); 
    assert(grid_.sizeg()>0);
    assert(grid_.inc(2)==1);
    assert(dis=='d' || dis=='g' || dis=='s');
    assert(src!=0);

    if( !grid_.active() )return;

    bc_[0]=px;
    bc_[1]=py;
    bc_[2]=pz;
    
    setup();

    const short nghosts=ghost_pt();

    int  incx_src=0, incy_src=0;
    if(dis=='g' || dis=='s'){
        incx_src=grid_.gdim(1)*grid_.gdim(2);
        incy_src=grid_.gdim(2);
    }else if(dis=='d'){
        incx_src=dim(2)*dim(1);
        incy_src=dim(2);
    }else{
        cout<<"Undefined distribution to build GridFunc<T>!!"<<endl;
        mype_env().globalExit(2);
    }

    int istart=0;
    int jstart=0;
    int kstart=0;
    if(dis=='g'){
        istart=mype_env().my_mpi(0)*dim(0);
        jstart=mype_env().my_mpi(1)*dim(1);
        kstart=mype_env().my_mpi(2)*dim(2);
    }else if(dis=='s'){ // "shifted" option
        istart=mype_env().my_mpi(0)*dim(0)+(grid_.gdim(0)>>1);
        jstart=mype_env().my_mpi(1)*dim(1)+(grid_.gdim(1)>>1);
        kstart=mype_env().my_mpi(2)*dim(2)+(grid_.gdim(2)>>1);
    }
 
    uu_ = new T[grid_.sizeg()];
    memset(uu_,0,grid_.sizeg()*sizeof(T));

    const size_t sdim2=dim_[2]*sizeof(T);

    if(dis=='g' || dis=='d')
    {
        for(int ix = 0;ix < dim_[0];ix++) {

            int ix1 = (ix+nghosts)*incx_+nghosts+nghosts*incy_;
            int ix2 = (ix+istart)*incx_src+jstart*incy_src;

            for(int iy = 0;iy < dim_[1];iy++) {

                int iy1 = ix1+iy*incy_;
                int iy2 = ix2+iy*incy_src;
 
                memcpy(uu_+iy1, src+kstart+iy2, sdim2);
            }
        }
    }
    else
    {
        const int gdim0=grid_.gdim(0);
        const int gdim1=grid_.gdim(1);
        const int gdim2=grid_.gdim(2);

        for(int ix = 0;ix < dim_[0];ix++) {

            int ix1 = (ix+nghosts)*incx_+nghosts;
            int ix2 = ( (ix+istart) % gdim0 )*incx_src ;

            for(int iy = 0;iy < dim_[1];iy++) {

                int iy1 = ix1+(iy+nghosts)*incy_;
                int iy2 = ix2+( (iy+jstart)% gdim1 )*incy_src;
 
                const T* __restrict__ psrc=src+iy2;
                T* __restrict__ puu=uu_+iy1;
                for(int iz = 0;iz < dim_[2];iz++) {
                    const int ii = (kstart+iz)%gdim2; 
                    puu[iz]=psrc[ii];
                }
            }

        }
    }
    
    updated_boundaries_=false;

    const short shift=ghost_pt();
    const int dimxy=(dim(1)+2*ghost_pt())*dim(0);
    
    // resize static buffers if needed
    if(mype_env().n_mpi_tasks()>1){
        int size_max = shift*grid_.inc(0);
        size_max = max(size_max, shift*dimxy);
        size_max = max(size_max, shift*dim(0)*grid_.inc(1));

        if( size_max>(int)buf1_.size() ){
            buf1_.resize(size_max);
            buf2_.resize(size_max);
            buf3_.resize(size_max);
            buf4_.resize(size_max);
        }
    }
}

template <typename T>
GridFunc<T>& GridFunc<T>::operator-=(const GridFunc<T>& func)
{
    assert( func.grid_.sizeg()==grid_.sizeg() );
    assert( func.grid_.ghost_pt()==grid_.ghost_pt() );
    assert( this!=&func );

    if( !grid_.active() )return *this;

    int n=grid_.sizeg();
    T minus=-1.;
    Taxpy(n, minus, func.uu_, uu_);
            
    updated_boundaries_=(func.updated_boundaries() && updated_boundaries_);

    return *this;    
}

template <typename T>
void GridFunc<T>::set_max(const T val)
{
    int n=grid_.sizeg();

    for (int i=0; i< n; i++){
        T alpha=fabs(uu_[i]);
        if( alpha>val )uu_[i]=(T)(val*uu_[i]/alpha);
    }

}


template <typename T>
void GridFunc<T>::setValues(const int n, const T* src, const int pos)
{
    assert((pos + n) <= grid_.sizeg());

    size_t ssize=n*sizeof(T);
    //int ione=1;
    memcpy(&uu_[pos], src, ssize);
}

template <typename T>
void GridFunc<T>::setValues(const GridFunc<T>& src)
{
    assert(src.grid_.sizeg()==grid_.sizeg());
    
    setValues(src.grid_.sizeg(), src.uu_);
}

template <typename T>
GridFunc<T>& GridFunc<T>::operator=(const GridFunc<T>& func)
{
    assert(func.grid_.sizeg()==grid_.sizeg());

    if(this==&func)return *this;

    bc_[0]=func.bc_[0];
    bc_[1]=func.bc_[1];
    bc_[2]=func.bc_[2];
    updated_boundaries_=func.updated_boundaries();

    setValues(func.grid_.sizeg(), func.uu_);
            
    return *this;    

}

template <typename T>
void GridFunc<T>::setValues(const T val)
{
    const int n=grid_.sizeg();
    T* __restrict__ pu=&uu_[0];

    for (int i=0; i< n; i++)
        pu[i] = val;
}

template <typename T>
GridFunc<T>& GridFunc<T>::operator=(const T val)
{
    setValues(val);

    updated_boundaries_=true;
    return *this;    
}

template <typename T>
GridFunc<T> GridFunc<T>::operator+(const GridFunc<T> &A)
{
    int     n = grid_.sizeg();

    assert(n==A.grid_.sizeg());
    GridFunc<T> tmp(A);
    T* __restrict__ tu=&tmp.uu_[0];
    const T* __restrict__ pu=&uu_[0];

    for (int i=0; i<n; i++)
        tu[i] += pu[i];
    tmp.updated_boundaries_=(A.updated_boundaries_ && updated_boundaries_);

    return tmp;
}

template <typename T>
GridFunc<T> GridFunc<T>::operator-(const GridFunc<T> &A)
{ 
    const int  n = grid_.sizeg();
    assert(n==A.grid_.sizeg());
    assert(grid_.size()==A.grid_.size());

    GridFunc<T> tmp(*this);
    const T* __restrict__ au=A.uu(0);
    T* __restrict__       tu=tmp.uu(0);

    for (int i=0; i<n; i++){
        tu[i] -= au[i];
    }

    tmp.updated_boundaries_=(A.updated_boundaries() && updated_boundaries_);

    return tmp;
}

template <typename T>
GridFunc<T> GridFunc<T>::operator*(const double val)
{
    int n=grid_.sizeg();    
    GridFunc<T> tmp(grid_,bc_[0],bc_[1],bc_[2]);

    T* __restrict__ tu=&tmp.uu_[0];
    const T* __restrict__ pu=&uu_[0];
    for (int i=0; i< n; i++)
        tu[i] = val*pu[i];

    tmp.set_updated_boundaries(updated_boundaries_);

    return tmp;    
}

template <typename T>
GridFunc<T>& GridFunc<T>::operator+=(T alpha)
{    
    const int n=grid_.sizeg();    
    T* __restrict__ pu=&uu_[0];
    for(int i=0;i<n;i++)
        pu[i]+= alpha;
            
    return *this;    
}

template <typename T>
GridFunc<T>& GridFunc<T>::operator-=(const T alpha)
{    
    int n=grid_.sizeg();    
    for(int i=0;i<n;i++)
        uu_[i]-= alpha;
            
    return *this;
}

template <typename T>
GridFunc<T>& GridFunc<T>::operator*=(const GridFunc<T>& B)
{    
    int n=grid_.sizeg();    
    assert(B.grid_.sizeg()==n);
    const T * __restrict__ v=B.uu_;
    T* __restrict__ pu=&uu_[0];
    for(int i=0;i<n;i++){
        pu[i]*=v[i];
    }
    updated_boundaries_=false;
                
    return *this;    
}

template <typename T>
GridFunc<T>& GridFunc<T>::operator/=(const GridFunc<T>& B)
{    
    int n=grid_.sizeg();    
    assert(B.grid_.sizeg()==n);
    T* const pu=&uu_[0];
    T* const bu=&B.uu_[0];

    for(int i=0;i<n;i++){
            pu[i]/=bu[i];
    }
    updated_boundaries_=false;
                
    return *this;    
}

template <typename T>
GridFunc<T>& GridFunc<T>::operator+=(const GridFunc<T>& func)
{
    assert(func.grid_.sizeg()==grid_.sizeg());
    
    if( !grid_.active() )return *this;

    Taxpy(func.grid_.sizeg(), 1., func.uu_, uu_);
                
    updated_boundaries_=(func.updated_boundaries() && updated_boundaries_);

    return *this;
}

template <typename T>
GridFunc<T>& GridFunc<T>::operator*=(const double alpha)
{    
    if( grid_.active() )
        MPscal(grid_.sizeg(), alpha, uu_);
                
    return *this;
}

template <typename T>
void GridFunc<T>::scal(const double alpha)
{
    if( !grid_.active() )return;
    
    MPscal(grid_.sizeg(), alpha, uu_);
}

template <typename T>
void GridFunc<T>::axpy(const double alpha, const GridFunc<T>& vv)
{
    assert(vv.grid_.sizeg()==grid_.sizeg());

    if( !grid_.active() )return;

    int n=grid_.sizeg();
    MPaxpy(n, alpha, vv.uu_, uu_);
                    
    updated_boundaries_=(vv.updated_boundaries() && updated_boundaries_);
}

template <typename T>
void GridFunc<T>::prod(const GridFunc<T>& A, const GridFunc<T>& B)
{
    assert(A.grid_.sizeg()==grid_.sizeg());
    assert(B.grid_.sizeg()==grid_.sizeg());

    prod_tm_.start();

    if( !grid_.active() )return;

    const T * const v1=A.uu_;
    const T * const v2=B.uu_;
    T* const pu=&uu_[0];
    const int n=grid_.sizeg();
    for(int i=0;i<n;i++)
    {
        pu[i]=(v1[i]*v2[i]);
    }

    updated_boundaries_=(A.updated_boundaries() && B.updated_boundaries());

    prod_tm_.stop();                    
}

template <typename T>
void GridFunc<T>::diff(const GridFunc<T>& A, const GridFunc<T>& B)
{
    assert(A.grid_.sizeg()==grid_.sizeg());
    assert(B.grid_.sizeg()==grid_.sizeg());

    if( !grid_.active() )return;

    const T * const v1=A.uu_;
    const T * const v2=B.uu_;
    T* const pu=&uu_[0];
    const int n=grid_.sizeg();
    for(int i=0;i<n;i++){
        pu[i]=(v1[i]-v2[i]);
    }

    updated_boundaries_=(A.updated_boundaries() && B.updated_boundaries());
}

// Initialize data with ghosts uu_ by copying data without ghosts vv
template <typename T>
void GridFunc<T>::assign(const double* const vv, const char dis)
{
    if( !grid_.active() )return;

    const short nghosts=ghost_pt();

    assert(grid_.inc(2)==1);
 
    if(vv!=0)
    {
        if( dis=='g' )
        {
            const int ldz=grid_.gdim(2);
            const int istart=mype_env().my_mpi(0)*dim_[0];
            const int jstart=mype_env().my_mpi(1)*dim_[1];
            const int kstart=mype_env().my_mpi(2)*dim_[2];
            const int gincx=grid_.gdim(1)*ldz;
            const int gincy=ldz;
        
            for(int ix = 0;ix < dim_[0];ix++)
            {
                int ix1 = (ix+ nghosts)*incx_;
                int ix2 = (ix+istart)*gincx;

                for(int iy = 0;iy < dim_[1];iy++)
                {
                    int iy1 = ix1+(iy+nghosts)*incy_;
                    int iy2 = ix2+(iy+jstart)*gincy+kstart;
                        
                    MPcpy(uu_+iy1+nghosts, vv+iy2, dim_[2]);
                } 
            }

        }
        else
        {
            const int incx2=dim_[2]*dim_[1];
            const int incy2=dim_[2];
            for(int ix = 0;ix < dim_[0];ix++)
            {
                int ix1 = (ix+ nghosts)*incx_;
                int ix2 = ix*incx2;

                const double* const pv=vv+ix2;
                T* pu=uu_+ix1+nghosts+nghosts*incy_;
                for(int iy = 0;iy < dim_[1];iy++)
                {
                    int iy1 = iy*incy_;
                    int iy2 = iy*incy2;
                        
                    MPcpy(pu+iy1, pv+iy2, dim_[2]);  
                } 
            }
        }

        updated_boundaries_=false;
    }
    else
    {
        updated_boundaries_=false;
    }

}

// Initialize data with ghosts uu_ by copying data without ghosts vv
template <typename T>
void GridFunc<T>::assign(const float* const vv, const char dis)
{
    if( !grid_.active() )return;

    const short nghosts=ghost_pt();

    if(vv!=0)
    {
        if( dis=='g' )
        {
            const int ldz=grid_.gdim(2);
            const int istart=mype_env().my_mpi(0)*dim_[0];
            const int jstart=mype_env().my_mpi(1)*dim_[1];
            const int kstart=mype_env().my_mpi(2)*dim_[2];
            const int gincx=grid_.gdim(1)*ldz;
            const int gincy=ldz;
        
            for(int ix = 0;ix < dim_[0];ix++)
            {
                int ix1 = (ix+ nghosts)*incx_;
                int ix2 = (ix+istart)*gincx;

                for(int iy = 0;iy < dim_[1];iy++)
                {
                    int iy1 = ix1+(iy+nghosts)*incy_;
                    int iy2 = ix2+(iy+jstart)*gincy+kstart;
                        
                    MPcpy(uu_+iy1+nghosts, vv+iy2, dim_[2]);
                } 
            }

        }
        else
        {
            const int incx2=dim_[2]*dim_[1];
            const int incy2=dim_[2];
            for(int ix = 0;ix < dim_[0];ix++)
            {
                int ix1 = (ix+ nghosts)*incx_;
                int ix2 = ix*incx2;

                const float* const pv=vv+ix2;
                T* pu=uu_+ix1+nghosts+nghosts*incy_;
                for(int iy = 0;iy < dim_[1];iy++)
                {
                    int iy1 = iy*incy_;
                    int iy2 = iy*incy2;
                        
                    MPcpy(pu+iy1, pv+iy2, dim_[2]);  
                } 
            }
        }

        updated_boundaries_=false;
    }
    else
    {
        updated_boundaries_=false;
    }

}

template <typename T>
int GridFunc<T>::count_threshold(const T threshold)
{
    if( !grid_.active() )return 0;

    const short nghosts=ghost_pt();

    trade_boundaries();
 
    int icount=0;

    for(int ix = 0;ix < dim_[0];ix++) {

        int ix1 = (ix+ nghosts)*incx_+nghosts;

        for(int iy = 0;iy < dim_[1];iy++) {

            int iy1 = ix1+(iy+ nghosts)*incy_;
                
            T* const pu=&uu_[iy1];
            for(int iz =0;iz < dim_[2];iz++) {
    
                if(fabs(pu[iz])>threshold)icount++;

            } 
        } 

    } 
    assert(icount<=size());
    
#ifdef USE_MPI
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if( mype_env().n_mpi_tasks()>1 ){
        int  sum=0;
        int rc = mmpi.allreduce(&icount, &sum, 1, MPI_SUM);
        if(rc!=MPI_SUCCESS){
            cout<<"MPI_Allreduce double sum failed!!!"<<endl;
            mype_env().globalExit(2);
        }
        icount=sum;
    }
#endif

    return icount;

}

template <typename T>
void GridFunc<T>::sqrt_func()
{
    if( !grid_.active() )return;

    const short shift=ghost_pt();

    assert(grid_.inc(2)==1);
 
    for(int ix = 0;ix < dim_[0];ix++) {

        int ix1 = (ix+ shift)*incx_+shift;

        for(int iy = 0;iy < dim_[1];iy++) {

            int iy1 = ix1+(iy+ shift)*incy_;    
            T* __restrict__ pu=&uu_[iy1];
        
            for(int iz =0;iz < dim_[2];iz++) {
    
                assert(pu[iz]>=0.);
                pu[iz]= sqrt(pu[iz]);

            } 
   
        } 

    } 

    updated_boundaries_=false;
}

template <typename T>
void GridFunc<T>::inv_sqrt()
{
    if( !grid_.active() )return;

    const short shift=ghost_pt();

    assert(grid_.inc(2)==1);
 
    for(int ix = 0;ix < dim_[0];ix++) {

        int ix1 = (ix+ shift)*incx_+shift;

        for(int iy = 0;iy < dim_[1];iy++) {

            int iy1 = ix1+(iy+ shift)*incy_;
            T* __restrict__ pu=&uu_[iy1];
            
            for(int iz =0;iz < dim_[2];iz++) {
    
                assert(pu[iz]>1.e-5);
                pu[iz]= 1./sqrt(pu[iz]);

            } 
   
        } 

    } 

    updated_boundaries_=false;
}

template <typename T>
void GridFunc<T>::inv()
{
    if( !grid_.active() )return;

    const short shift=ghost_pt();

    assert(grid_.inc(2)==1);
 
    for(int ix = 0;ix < dim_[0];ix++) {

        int ix1 = (ix+ shift)*incx_+shift;

        for(int iy = 0;iy < dim_[1];iy++) {

            int iy1 = ix1+(iy+ shift)*incy_;
            T* __restrict__ pu=&uu_[iy1];
            
            for(int iz =0;iz < dim_[2];iz++) {
    
                assert(pu[iz]>1.e-8);
                pu[iz]= 1./pu[iz];

            } 
   
        } 

    } 

    updated_boundaries_=false;
}

// Print the values of the function in 2 columns
// 1st column: value to use for the ordering
// 2nd column: value of the function
// Transfer datas to task 0 and print from task 0
template <typename T>
void GridFunc<T>::print(ostream& tfile)
{
    if( !grid_.active() )return;

    const short nghosts=ghost_pt();

    assert(grid_.inc(2)==1);

    tfile<<mype_env().n_mpi_task(0)*dim_[0]<<"\t"
         <<mype_env().n_mpi_task(1)*dim_[1]<<"\t"
         <<mype_env().n_mpi_task(2)*dim_[2]<<endl;
    tfile<<grid_.hgrid(0)<<"\t"<<grid_.hgrid(1)<<"\t"<<grid_.hgrid(2)<<endl;

    double  size_x =grid_.hgrid(1)*grid_.hgrid(2)*dim_[1]*dim_[2];
    double  size_y= grid_.hgrid(0)*grid_.hgrid(2)*dim_[2];
    double  size_z= grid_.hgrid(0)*grid_.hgrid(1);

    T* work=new T[grid_.sizeg()];
    double  start[3]={grid_.start(0),grid_.start(1),grid_.start(2)};
    
    T* vv=uu_;

    if(uu_!=0){
#ifdef USE_MPI
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        int     mpirc;
        
        for(int i=0;i<mype_env().n_mpi_tasks();i++)
        {
            if(i>0)
            {
            if(i==mytask_){
                mpirc=mmpi.send(uu_, grid_.sizeg(), 0);
                mpirc=mmpi.send(start, 3, 0);
                cout<<" Print GridFunc<T> from task "<<mytask_<<endl;
            }else{
                if(0==mytask_){
                    mpirc=mmpi.recv(work, grid_.sizeg(), i);
                    if(mpirc!=MPI_SUCCESS)
                        cout<<"print: MPI_Recv work failed!!!"<<endl;
                    mpirc=mmpi.recv(start, 3, i);                
                    if(mpirc!=MPI_SUCCESS)
                        cout<<"print: MPI_Recv start[] failed!!!"<<endl;
                    vv=work;
                }
            }
            }
#endif
            if(mype_env().onpe0()){
                for(int ix = 0;ix < dim_[0];ix++) {

                    int ix1 = (ix+ nghosts)*incx_+nghosts;

                    for(int iy = 0;iy < dim_[1];iy++) {

                        int iy1 = ix1+(iy+ nghosts)*incy_;
 
                        for(int iz =0;iz < dim_[2];iz++) {
 
                            int iz1 = iy1+iz;

                            tfile<<(start[0]+ix*grid_.hgrid(0))*size_x
                              +(start[1]+iy*grid_.hgrid(1))*size_y
                              +(start[2]+iz*grid_.hgrid(2))*size_z<<"\t"
                             <<vv[iz1]<<endl;

                        }
 
                    }

                }
            }

#ifdef USE_MPI
        }
#endif
    }
    delete[] work;


}

// Print the values of the function in formatted .plt file for gOpenMol
template <typename T>
void GridFunc<T>::write_plt(const char str[])const
{
    if( !grid_.active() )return;

    ofstream  tfile;
    tfile.open(str);

    tfile<<3<<" "<<2<<endl;

    tfile<<dim(0)<<" "<<dim(1)<<" "<<dim(2)<<endl;

    tfile<<0.<<" "<<(dim(0)-1)*grid_.hgrid(0)<<" ";
    tfile<<0.<<" "<<(dim(1)-1)*grid_.hgrid(1)<<" ";
    tfile<<0.<<" "<<(dim(2)-1)*grid_.hgrid(2)<<endl;

    write_zyx(tfile);

}

template <typename T>
void GridFunc<T>::write_xyz(ofstream&  tfile)const
{
    if( !grid_.active() )return;

    const short shift=ghost_pt();

    assert(grid_.inc(2)==1);


    const T * const vv=uu_;

    if(uu_!=0){
        for(int ix = 0;ix < dim(0);ix++) {
 
            int ix1 = (ix+shift)*incx_+ shift;

            for(int iy = 0;iy < dim(1);iy++) {

                int iy1 = ix1+(iy+ shift)*incy_;
 
                for(int iz =0;iz < dim(2);iz++){

                    int iz1 = iy1+iz;

                    tfile<<vv[iz1]<<endl;

                }
 
            }

        }

    }

}

template <typename T>
void GridFunc<T>::global_xyz_task0(T *global_func)
{
    if( !grid_.active() )return;

    const short shift=ghost_pt();

    T* work=new T[grid_.sizeg()];

    const int gincx=mype_env().n_mpi_task(1)*dim(1)
                   *mype_env().n_mpi_task(2)*dim(2);
    const int gincy=mype_env().n_mpi_task(2)*dim(2);


    for(int ii=0;ii<dim(0);ii++)
    for(int jj=0;jj<dim(1);jj++)
            memcpy(global_func+ii*gincx+jj*gincy,
               uu_+(ii+shift)*incx_+(jj+shift)*incy_+shift,
               dim(2)*sizeof(T));

#ifdef USE_MPI
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int mpirc;

    for(int i=0;i<mype_env().n_mpi_task(0);i++)
    for(int j=0;j<mype_env().n_mpi_task(1);j++)
    for(int k=0;k<mype_env().n_mpi_task(2);k++){
        int istart=i*dim(0);
            int jstart=j*dim(1);
        int kstart=k*dim(2);
        int ipe=mype_env().xyz2task(i,j,k);

        if(ipe>0)
        {
        if(mype_env().onpe0()){
            mpirc=mmpi.recv(work, grid_.sizeg(), ipe);
            if(mpirc!=MPI_SUCCESS)
                    cout<<"print: MPI_Recv work failed!!!"<<endl;
            for(int ii=0;ii<dim(0);ii++)
            for(int jj=0;jj<dim(1);jj++)
                 memcpy(global_func+(istart+ii)*gincx+(jstart+jj)*gincy+kstart,
                        work+(ii+shift)*incx_+(jj+shift)*incy_+shift,
                        dim(2)*sizeof(T));
        }else if(mytask_==ipe){
            mpirc=mmpi.send(uu_, grid_.sizeg(), 0);
        }
        }
        
        MPI_Barrier( mype_env().comm() );        
    }
#endif
    
    delete[] work;

}

template <typename T>
void GridFunc<T>::write_global_xyz(ofstream&  tfile)
{
    if( !grid_.active() )return;

    assert(grid_.inc(2)==1);

    const short shift=ghost_pt();

    const long gsize=grid_.gsize();
    T* global_func=new T[gsize];
    T* work=new T[grid_.sizeg()];

    const int gincx=mype_env().n_mpi_task(1)*dim(1)
                   *mype_env().n_mpi_task(2)*dim(2);
    const int gincy=mype_env().n_mpi_task(2)*dim(2);

    size_t sizez=dim(2)*sizeof(T);
    const int ldim0=dim(0);
    const int ldim1=dim(1);

    for(int ii=0;ii<ldim0;ii++)
    for(int jj=0;jj<ldim1;jj++)
            memcpy(global_func+ii*gincx+jj*gincy,
               uu_+(ii+shift)*incx_+(jj+shift)*incy_+shift,
               sizez);

#ifdef USE_MPI
    if(mype_env().onpe0())
            cout<<"GridFunc<T>::write_global_xyz, Collect data on PE 0"<<endl;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int mpirc;
    const int ldim2=dim(2);

    for(int i=0;i<mype_env().n_mpi_task(0);i++)
    for(int j=0;j<mype_env().n_mpi_task(1);j++)
    for(int k=0;k<mype_env().n_mpi_task(2);k++){
        const int istart=i*ldim0;
            const int jstart=j*ldim1;
        const int kstart=k*ldim2;
        const int ipe=mype_env().xyz2task(i,j,k);

        if(ipe>0)
        {
        if(mytask_==0)
        {
            mpirc=mmpi.recv(work, grid_.sizeg(), ipe);
            if(mpirc!=MPI_SUCCESS)
                cout<<"GridFunc<T>::write_global_xyz, MPI_Recv work failed!!!"<<endl;
            for(int ii=0;ii<ldim0;ii++)
            for(int jj=0;jj<ldim1;jj++)
                 memcpy(global_func+(istart+ii)*gincx+(jstart+jj)*gincy+kstart,
                        work+(ii+shift)*incx_+(jj+shift)*incy_+shift,
                        sizez);
        }else if(mytask_==ipe){
            mpirc=mmpi.send(uu_, grid_.sizeg(), 0);
            if(mpirc!=MPI_SUCCESS)
                    cout<<"GridFunc<T>::write_global_xyz, MPI_Send work failed!!!"<<endl;
        }
        }
        
        MPI_Barrier( mype_env().comm() );        
    }
#endif
    
    delete[] work;

    // write from PE 0 only
    if(mype_env().onpe0()){
            cout<<"GridFunc<T>::write_global_xyz, Write data from PE 0"<<endl;
        for(long ix = 0;ix < gsize;ix++){
                tfile<<global_func[ix]<<endl;
            }
    }
#ifdef USE_MPI
        MPI_Barrier( mype_env().comm() );        
#endif
    
    delete[] global_func;

}

template <typename T>
void GridFunc<T>::write_global_x(const char str[])
{
    if( !grid_.active() )return;

    ofstream  tfile(str);

    T* global_func=new T[grid_.gsize()];

    this->global_xyz_task0(global_func);

    if(mype_env().onpe0()){
        int incx=grid_.gdim(1)*grid_.gdim(2);
        for(int ix = 0;ix < grid_.gdim(0);ix++){
            T alpha=0.;
            for(int iy = 0;iy < grid_.gdim(1);iy++)
            for(int iz = 0;iz < grid_.gdim(2);iz++){
                alpha+=global_func[ix*incx
                                  +iy*grid_.gdim(2)+iz];
            }
            tfile<<alpha/((T)incx_)<<endl;
        }

    }
    
    delete[] global_func;

}

#ifdef USE_MPI
template <typename T>
void GridFunc<T>::allGather(T *global_func)const
{
    assert(grid_.inc(2)==1);

    all_gather_tm_.start();
    
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    
    //if(mype_env().onpe0())
    //    cout<<"GridFunc<T>::allGather()"<<endl;
    if( !grid_.active() )return;
    
    const int gincx=grid_.gdim(1)*grid_.gdim(2);
    const int gincy=grid_.gdim(2);

    const short shift=ghost_pt();

    const int ldim[3]={dim(0),dim(1),dim(2)};

    if( mype_env().n_mpi_tasks()==1 ){
        for(int ii=0;ii<ldim[0];ii++)
        for(int jj=0;jj<ldim[1];jj++)
             memcpy(global_func+ii*gincx+jj*gincy,
                        uu_+(ii+shift)*incx_+(jj+shift)*incy_+shift,
                        ldim[2]*sizeof(T));
        return;
    }
    
    // Compute and communicate displacements (used in MPI_Allgather)

    const int ntasks=mype_env().n_mpi_tasks();
    int* displs=new int[ntasks];
    int mpi_err;
#if 0 // old version
    const int istart=grid_.istart(0);
    const int jstart=grid_.istart(1);
    const int kstart=grid_.istart(2);
    int mydispl=istart*gincx+jstart*gincy+kstart;
    //if(mype_env().onpe0())
    //    cout<<"MPI_Allgather()"<<endl;
    mpi_err=MPI_Allgather(&mydispl, 1, MPI_INT,
                          displs,   1, MPI_INT, mype_env().comm());
    if(mpi_err!=MPI_SUCCESS)
        cout<<"GridFunc<T>::gather --- MPI_Allgather() call failed"<<endl;
#else
    for(int i=0;i<ntasks;i++){
        int other_start[3];
        for(short dir=0;dir<3;dir++){
            other_start[dir]=mype_env().other_tasks_dir(dir,i)*ldim[dir];
        }
        displs[i]=other_start[0]*gincx
                 +other_start[1]*gincy
                 +other_start[2];
        assert( displs[i]<grid_.gsize() );
    }
#endif

    const int sizeg=grid_.sizeg();    
    const int gsize = sizeg*ntasks;
    const size_t sdim2=ldim[2]*sizeof(T);
#if 0 // old version
    T* buffer=new T[sizeg];
    T* tbuffer;
    for(int i=0;i<ntasks;i++){
        if(mytask_==i){
            tbuffer=uu_;
        }else{
            tbuffer=buffer;
        }
        mmpi.bcast(tbuffer, sizeg, i);
        for(int ii=0;ii<ldim[0];ii++){
            T* dest=global_func+displs[i]+ii*gincx;
            T* src=tbuffer+(ii+shift)*incx_+shift;
            for(int jj=0;jj<ldim[1];jj++){
                memcpy(dest+jj*gincy,
                        src+(jj+shift)*incy,
                        sdim2);
            }
        }
    }
#else
    T* buffer=new T[gsize];    
    mpi_err=mmpi.allGather(uu_, sizeg, buffer, gsize);

    for(int i=0;i<ntasks;i++)
    for(int ii=0;ii<ldim[0];ii++){
        T* const dest=global_func+displs[i]+ii*gincx;
        const T* const src=buffer+i*sizeg+(ii+shift)*incx_+shift*(1+incy_);
        for(int jj=0;jj<ldim[1];jj++){
            memcpy(dest+jj*gincy,
                   src+jj*incy_, sdim2);
        }
    }
#endif
    delete[] buffer;
    delete[] displs;

    all_gather_tm_.stop();
}
#endif

// gather GridFunc<T> into T* on PE 0
template <typename T>
void GridFunc<T>::gather(T *global_func)const
{
    assert(grid_.inc(2)==1);

    gather_tm_.start();
    
    //if(mype_env().onpe0())
    //    cout<<"GridFunc<T>::gather()"<<endl;
    if( !grid_.active() )return;
    
    const int gincx=grid_.gdim(1)*grid_.gdim(2);
    const int gincy=grid_.gdim(2);

    const short shift=ghost_pt();

    const int ldim[3]={dim(0),dim(1),dim(2)};

    if( mype_env().n_mpi_tasks()==1 ){
        for(int ii=0;ii<ldim[0];ii++)
        for(int jj=0;jj<ldim[1];jj++)
             memcpy(global_func+ii*gincx+jj*gincy,
                        uu_+(ii+shift)*incx_+(jj+shift)*incy_+shift,
                        ldim[2]*sizeof(T));
        return;
    }
    
#ifdef USE_MPI
    // Compute and communicate displacements (used in MPI_Allgather)
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int* displs=NULL;
    const int ntasks=mype_env().n_mpi_tasks();
    if( mype_env().onpe0() ){
        displs=new int[ntasks];

        for(int i=0;i<ntasks;i++){
            int other_start[3];
            for(short dir=0;dir<3;dir++){
                other_start[dir]=mype_env().other_tasks_dir(dir,i)*ldim[dir];
            }
            displs[i]=other_start[0]*gincx
                     +other_start[1]*gincy
                     +other_start[2];
            assert( displs[i]<grid_.gsize() );
        }
    }

    const int sizeg=grid_.sizeg();
    const int gsize=sizeg*ntasks;

    T* buffer=NULL;
    if( mype_env().onpe0() )
        buffer=new T[gsize];
    
    mmpi.gather(uu_, sizeg, buffer, gsize, 0);

    if( mype_env().onpe0() ){
        const size_t sdim2=ldim[2]*sizeof(T);
        for(int i=0;i<ntasks;i++){
            T* const pdst = global_func+displs[i];
            T* const psrc = buffer+i*sizeg;
            for(int ii=0;ii<ldim[0];ii++){
                T* const       dst=pdst+ii*gincx;
                const T* const src=psrc+(ii+shift)*incx_+shift*(1+incy_);
                for(int jj=0;jj<ldim[1];jj++){
                    memcpy(dst+jj*gincy,
                           src+jj*incy_, sdim2);
                }
            }
        }

        delete[] buffer;
        delete[] displs;
    }
#endif

    gather_tm_.stop();
}

// scatter GridFunc<T> from PE 0 to all PEs
template <typename T>
void GridFunc<T>::scatterFrom(const GridFunc<T>& src)
{
    assert(grid_.inc(2)==1);

    scatter_tm_.start();
    
    //if(mype_env().onpe0())
    //    cout<<"GridFunc<T>::scatterFrom()"<<endl;
    if( !grid_.active() )return;
    
    const int incx_src=src.grid_.inc(0);
    const int incy_src=src.grid_.inc(1);

    const int ghosts_src=src.ghost_pt();

    const int ghosts_dst=ghost_pt();

    const int ldim[3]={dim(0),dim(1),dim(2)};

    if( mype_env().n_mpi_tasks()==1 ){
        for(int ii=0;ii<ldim[0];ii++)
        for(int jj=0;jj<ldim[1];jj++)
             memcpy(uu_+(ii+ghosts_dst)*incx_+(jj+ghosts_dst)*incy_+ghosts_dst,
                        src.uu_+(ii+ghosts_src)*incx_src+(jj+ghosts_src)*incy_src+ghosts_src,
                        ldim[2]*sizeof(T));
        return;
    }
    
#ifdef USE_MPI
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    // Compute and communicate displacements (used in MPI_Scatter)
    int* displs=NULL;
    const int ntasks=mype_env().n_mpi_tasks();
    if( mype_env().onpe0() ){
        displs=new int[ntasks];

        for(int i=0;i<ntasks;i++){
            int other_start[3];
            for(short dir=0;dir<3;dir++){
                other_start[dir]=mype_env().other_tasks_dir(dir,i)*ldim[dir];
            }
            displs[i]=other_start[0]*incx_src
                     +other_start[1]*incy_src
                     +other_start[2];
            assert( displs[i]<grid_.gsize() );
        }
    }

    const int sizeg=grid_.sizeg();
    const int gsize=sizeg*ntasks;

    T* buffer=NULL;
    if( mype_env().onpe0() )
        buffer=new T[gsize];
    
    if( mype_env().onpe0() ){ // fill buffer
        const size_t sdim2=ldim[2]*sizeof(T);
        for(int i=0;i<ntasks;i++){
            T* const       pdst = buffer+i*sizeg;
            const T* const psrc = src.uu_+displs[i];
            for(int ii=0;ii<ldim[0];ii++){
                T* const       ldst=pdst+(ii+ghosts_dst)*incx_   +ghosts_dst*(1+incy_);
                const T* const lsrc=psrc+(ii+ghosts_src)*incx_src+ghosts_src*(1+incy_src);
                for(int jj=0;jj<ldim[1];jj++){
                    memcpy(ldst+jj*incy_,
                           lsrc+jj*incy_src, sdim2);
                }
            }
        }
    }
    
    mmpi.scatter(buffer, gsize, uu_, sizeg, 0);

    if( mype_env().onpe0() ){
        delete[] buffer;
        delete[] displs;
    }
#endif   

    updated_boundaries_=false;

    scatter_tm_.stop();
}

// Initialize a global array from a GridFunc<T> object shifted by half
// the global grid size
template <typename T>
void GridFunc<T>::init_vect_shift(T *global_func)const
{
    if( !grid_.active() )return;

    const short shift=ghost_pt();

    assert(grid_.inc(2)==1);

    const int gincx=grid_.gdim(1)*grid_.gdim(2);
    const int gincy=grid_.gdim(2);


    // Compute and communicate displacements (used in MPI_Allgatherv)
    const int il=(grid_.gdim(0)>>1);
    const int jl=(grid_.gdim(1)>>1);
    const int kl=(grid_.gdim(2)>>1);
    const int dimil2=(grid_.dim(0)>>1);
    const int dimjl2=(grid_.dim(1)>>1);
    const int dimkl2=(grid_.dim(2)>>1);
    
    MGmol_MPI& mmpi=*(MGmol_MPI::instance());
    const int bufsize = mype_env().n_mpi_tasks();
    
    vector< vector<int> >  displx;
    displx.resize(2);
    displx[0].resize(bufsize);
    displx[1].resize(bufsize);
    
    vector< vector<int> >  disply;
    disply.resize(2);
    disply[0].resize(bufsize);
    disply[1].resize(bufsize);
    
    vector< vector<int> >  displz;
    displz.resize(2);
    displz[0].resize(bufsize);
    displz[1].resize(bufsize);
    
    int istart=(mype_env().my_mpi(0)*dim(0)+il)%grid_.gdim(0);
    int jstart=(mype_env().my_mpi(1)*dim(1)+jl)%grid_.gdim(1);
    int kstart=(mype_env().my_mpi(2)*dim(2)+kl)%grid_.gdim(2);
    mmpi.allGather(&istart, 1, &displx[0][0],bufsize);
    mmpi.allGather(&jstart, 1, &disply[0][0],bufsize);
    mmpi.allGather(&kstart, 1, &displz[0][0],bufsize);        

    istart=(mype_env().my_mpi(0)*dim(0)+dimil2+il)%grid_.gdim(0);
    jstart=(mype_env().my_mpi(1)*dim(1)+dimjl2+jl)%grid_.gdim(1);
    kstart=(mype_env().my_mpi(2)*dim(2)+dimkl2+kl)%grid_.gdim(2);
    mmpi.allGather(&istart, 1, &displx[1][0],bufsize);
    mmpi.allGather(&jstart, 1, &disply[1][0],bufsize);
    mmpi.allGather(&kstart, 1, &displz[1][0],bufsize); 

    int sizeg=grid_.sizeg();

    T* buffer=new T[sizeg];
    T* tbuffer;
    for(int i=0;i<mype_env().n_mpi_tasks();i++){
        if(mytask_==i){
            tbuffer=uu_;
        }else{
            tbuffer=buffer;
        }
        mmpi.bcast(tbuffer, sizeg, i);
        for(int a=0;a<2;a++)
        for(int b=0;b<2;b++){
            int displ1=displx[a][i]*gincx+disply[b][i]*gincy+displz[0][i];
            int displ2=displx[a][i]*gincx+disply[b][i]*gincy+displz[1][i];
            for(int ii=a*dimil2;ii<(a+1)*dimil2;ii++)
            for(int jj=b*dimjl2;jj<(b+1)*dimjl2;jj++){
                    memcpy(global_func+displ1+(ii-a*dimil2)*gincx
                                         +(jj-b*dimjl2)*gincy,
                           tbuffer+(ii+shift)*incx_+(jj+shift)*incy_+shift,
                           dimkl2*sizeof(T));
                    memcpy(global_func+displ2+(ii-a*dimil2)*gincx
                                         +(jj-b*dimjl2)*gincy,
                           tbuffer+(ii+shift)*incx_+(jj+shift)*incy_+shift+dimkl2,
                           dimkl2*sizeof(T));
            }
        }
    }
    delete[] buffer;
}

template <typename T>
void GridFunc<T>::write_zyx(ofstream&  tfile)const
{
    if( !grid_.active() )return;

    const short shift=ghost_pt();

    const T * const vv=uu_;

    if(uu_!=0){
        for(int iz =0;iz < dim(2);iz++){

            int iz1 = (iz+ shift);

            for(int iy = 0;iy < dim(1);iy++) {

                int iy1 = iz1+(iy+ shift)*incy_;
 
                for(int ix = 0;ix < dim(0);ix++) {
 
                    int ix1 = iy1+(ix+shift)*incx_;

                    tfile<<vv[ix1]<<endl;

                }
 
            }

        }

    }


}

// assign vv with values in member uu_
template <typename T>
void GridFunc<T>::init_vect(T* vv, const char dis)const
{
    if( !grid_.active() )return;

    assert(grid_.inc(2)==1);
    assert(vv!=0);
    
#ifdef USE_MPI
    if(dis=='g'){
        allGather(vv);
    }else
#else
    (void) dis; // unused
#endif    
    {
        getValues(vv);
    }
    //cout<<"Copy constructor for function on grid "<<grid_.level()<<endl;    

}

// assign vv with values in member uu_
template <typename T>
void GridFunc<T>::getValues(float* vv)const
{
    if( !grid_.active() )return;

    assert(grid_.inc(2)==1);
    assert(vv!=0);
    
    const short nghosts=ghost_pt();

    const int incx_dest=dim_[2]*dim_[1];
    const int incy_dest=dim_[2];
    const int ix0=nghosts*(incx_+1);

    for(int ix = 0;ix < dim_[0];ix++)
    {
        int ix1 = ix0 + ix*incx_;
        int ix2 = ix*incx_dest;

        float* pdest=vv+ix2;
        const T* pu=uu_+ix1+nghosts*incy_;
        
        for(int iy = 0;iy < dim_[1];iy++)
        {
            MPcpy(pdest, pu, dim_[2]);
            
            pdest += incy_dest;
            pu += incy_;
   
        } 

    } 
}

// assign vv with values in member uu_
template <typename T>
void GridFunc<T>::getValues(double* vv)const
{
    if( !grid_.active() )return;

    assert(grid_.inc(2)==1);
    assert(vv!=0);
  
    const short nghosts=ghost_pt();

    const int incx_dest=dim_[2]*dim_[1];
    const int incy_dest=dim_[2];
    const int ix0=nghosts*(incx_+1);

    for(int ix = 0;ix < dim_[0];ix++)
    {
        int ix1 = ix0 + ix*incx_;
        int ix2 = ix*incx_dest;

        double* pdest=vv+ix2;
        const T* pu=uu_+ix1+nghosts*incy_;
        
        for(int iy = 0;iy < dim_[1];iy++)
        {
            MPcpy(pdest, pu, dim_[2]);
            
            pdest += incy_dest;
            pu += incy_;
   
        } 

    } 
}

template <typename T>
void GridFunc<T>::initTrigo3d(const short bc[3],
                           const int n[3])
{
    if( !grid_.active() )return;
    
    const short nghosts = ghost_pt();

    const int ilow = mype_env().my_mpi(0) * dim_[0];
    const int jlow = mype_env().my_mpi(1) * dim_[1];
    const int klow = mype_env().my_mpi(2) * dim_[2];

    const bool last[3]={(mype_env().my_mpi(0)==(mype_env().n_mpi_task(0)-1)),
                        (mype_env().my_mpi(1)==(mype_env().n_mpi_task(1)-1)),
                        (mype_env().my_mpi(2)==(mype_env().n_mpi_task(2)-1))};

    double (*f0)(double);
    double (*f1)(double);
    double (*f2)(double);
    int init[3]={0,0,0};
    int end[3]={dim_[0]+2*nghosts,dim_[1]+2*nghosts,dim_[2]+2*nghosts};
    if (bc[0]==1){
        f0 = &cos;
    }else{
        f0 = &sin;
    }
    if (bc[0]==0 && ilow==0){
        init[0]=nghosts;
    }
    if (bc[0]==0 && last[0]){
        end[0]-=nghosts;
    }
    if (bc[1]==1){
        f1 = &cos;
    }else{
        f1 = &sin;
    }
    if (bc[1]==0 && jlow==0){
        init[1]=nghosts;
    }
    if (bc[1]==0 && last[1]){
        end[1]-=nghosts;
    }
    if (bc[2]==1){
        f2 = &cos;
    }else{
        f2 = &sin;
    }
    if (bc[2]==0 && klow==0){
        init[2]=nghosts;
    }
    if (bc[2]==0 && last[2]){
        end[2]-=nghosts;
    }
    
    const double h0 = grid_.hgrid(0);
    const double h1 = grid_.hgrid(1);
    const double h2 = grid_.hgrid(2);
    
    const double k0 = (double)n[0]*2.*M_PI/grid_.ll(0);
    const double k1 = (double)n[1]*2.*M_PI/grid_.ll(1);
    const double k2 = (double)n[2]*2.*M_PI/grid_.ll(2);
    
    for(int ix = init[0];ix < end[0];ix++) {

        double xc=h0*(ix-nghosts+ilow);

        int iix = ix*grid_.inc(0);

        for(int iy = init[1];iy < end[1];iy++) {

            double yc=h1*(iy-nghosts+jlow);

            int iiy=iy*grid_.inc(1)+iix;

            for(int iz = init[2];iz < end[2];iz++) {
            
                double zc=h2*(iz-nghosts+klow);

                int iiz=iiy+iz;

                uu_[iiz]= (T)((*f0)(xc*k0)*(*f1)(yc*k1)*(*f2)(zc*k2));

            } 
       
        } 

    } 

    updated_boundaries_=false;    
}
 
template <typename T>
void GridFunc<T>::initCos3d(const T k[3])
{
    if( !grid_.active() )return;
    
    const short nghosts = ghost_pt();

    const int ilow = grid_.istart(0);
    const int jlow = grid_.istart(1);
    const int klow = grid_.istart(2);

    int init[3]={0,0,0};
    int end[3]={dim_[0],dim_[1],dim_[2]};

    const double h[3] = {grid_.hgrid(0),grid_.hgrid(1),grid_.hgrid(2)};
    
    for(int ix = init[0];ix < end[0];ix++) {

        const double xc=h[0]*(ix+ilow);

        const int iix = (ix+nghosts)*grid_.inc(0);

        for(int iy = init[1];iy < end[1];iy++) {

            const double yc=h[1]*(iy+jlow);

            const int iiy=(iy+nghosts)*grid_.inc(1)+iix;

            for(int iz = init[2];iz < end[2];iz++) {
            
                const double zc=h[2]*(iz+klow);

                const int iiz=iiy+iz+nghosts;

                uu_[iiz]= (T)(cos(xc*k[0])*cos(yc*k[1])*cos(zc*k[2]));

            } 
       
        }

    } 

    updated_boundaries_=false;    
}

template <typename T> 
void GridFunc<T>::initTriLin(const T a[4], const bool wghosts)
{
    //cout<<"a[0]="<<a[0]<<",a[1]="<<a[1]<<",a[2]="<<a[2]<<",a[3]="<<a[3]<<endl;
    if( !grid_.active() )return;
    
    const short nghosts = ghost_pt();

    const unsigned dim[3]={grid_.dim(0),grid_.dim(1),grid_.dim(2)};

    const int ilow = grid_.istart(0);
    const int jlow = grid_.istart(1);
    const int klow = grid_.istart(2);

    int init[3];
    int end[3];
    if( wghosts ){
        for(int i=0;i<3;i++){
            init[i]=0;
            end[i] =dim[i]+2*nghosts;
        }
    }
    else
        for(int i=0;i<3;i++){
            init[i]=nghosts;
            end[i] =dim[i]+nghosts;
        }
    //cout<<"end[0]="<<end[0]<<",end[1]="<<end[1]<<",end[2]="<<end[2]<<endl;
    
    const double h[3] = {grid_.hgrid(0),grid_.hgrid(1),grid_.hgrid(2)};
    
    for(int ix = init[0];ix < end[0];ix++) {

        const double xc=h[0]*(ix-nghosts+ilow)+grid_.origin(0);

        const int iix = ix*grid_.inc(0);

        for(int iy = init[1];iy < end[1];iy++) {

            const double yc=h[1]*(iy-nghosts+jlow)+grid_.origin(1);

            const int iiy=iy*grid_.inc(1)+iix;

            for(int iz = init[2];iz < end[2];iz++) {
            
                const double zc=h[2]*(iz-nghosts+klow)+grid_.origin(2);

                const int iiz=iiy+iz;

                uu_[iiz]= (T)(a[0]+xc*a[1]+yc*a[2]+zc*a[3]);
            } 
       
        }

    } 

    if( wghosts )updated_boundaries_=true;
    else         updated_boundaries_=false;
}

template <typename T> 
void GridFunc<T>::init_radial(const double center[3],
                              const double rcore,
                              const short ftype)
{
    if( !grid_.active() )return;
    
    if( mype_env().onpe0() )
        cout<<"Init_radial() at center ("<<center[0]<<","
                                         <<center[1]<<","
                                         <<center[2]<<")"<<endl;

    int     iix, iiy, iiz;
    
    T  r, r2;
    
    const double h0 = grid_.hgrid(0);
    const double h1 = grid_.hgrid(1);
    const double h2 = grid_.hgrid(2);
    
    const int ilow = mype_env().my_mpi(0) * grid_.dim(0);
    const int jlow = mype_env().my_mpi(1) * grid_.dim(1);
    const int klow = mype_env().my_mpi(2) * grid_.dim(2);

    const short nghosts = ghost_pt();

    for(int ix = 0;ix < dim_[0];ix++) {

        double xc=h0*(ix+ilow);

        iix = (ix+nghosts)*grid_.inc(0);

        for(int iy = 0;iy < dim_[1];iy++) {

            double yc=h1*(iy+jlow);

            iiy=(iy+nghosts)*grid_.inc(1)+iix;

            for(int iz = 0;iz < dim_[2];iz++) {
            
                double zc=h2*(iz+klow);

                iiz=iiy+iz+nghosts;

                double x = xc - center[0];
                double y = yc - center[1];
                double z = zc - center[2];

                r2 = (T)(x*x + y*y + z*z);
                r  = sqrt(r2);
                uu_[iiz]=radial_func(r, rcore, ftype);
                // uu_[iiz]= (*rfunc)(r,rcore);
                //uu_[iiz]= 1.;

            } 
       
        } 

    } 

    updated_boundaries_=false;    

    if( mype_env().onpe0() )
        cout<<" init_radial done"<<endl;
}

template <typename T> 
void GridFunc<T>::print_radial(const char str[])
{
    if( !grid_.active() )return;

    char   filename[30], extension[10];
    strcpy(filename,str);
    sprintf(extension, "%d", mype_env().my_mpi(2));
    strcat(filename,extension);

    ofstream  tfile(filename);
    
    const int gpt=ghost_pt();

    //if(mype_env().my_mpi(2)>0)return;
    if(mype_env().my_mpi(1)!=(mype_env().n_mpi_task(1)>>1))return;
    if(mype_env().my_mpi(0)!=(mype_env().n_mpi_task(0)>>1))return;
    
    //cout << " Print radial function...\n";
     
    int ix=gpt;
    int iy=gpt;
    
    if(mype_env().n_mpi_task(0)%2)ix+=(dim(0)>>1);
    if(mype_env().n_mpi_task(1)%2)iy+=(dim(1)>>1);
    
    int iix = ix*grid_.inc(0)+iy*grid_.inc(1);

    int izend=gpt+dim(2);
    if(gpt>0)izend++;
    for(int iz = gpt;iz < izend;iz++) {

        int iiz = iix+iz;

        tfile<<grid_.start(2)+(iz-gpt)*grid_.hgrid(2)
             <<"\t"<<uu_[iiz]<<endl;

    } 
} 

template <typename T>
void GridFunc<T>::init_rand()
{
    if( !grid_.active() )return;

#ifdef DEBUG
    if( mype_env().onpe0() )
        cout << " Begin init_rand ...\n";
#endif
    const int gpt=ghost_pt();
    
    const int xoff = grid_.istart(0);
    const int yoff = grid_.istart(1);
    const int zoff = grid_.istart(2);

    T* xrand=new T[grid_.gdim(0)];
    T* yrand=new T[grid_.gdim(1)];
    T* zrand=new T[grid_.gdim(2)];

    // Generate x, y, z random number sequences
    for(int idx = 0;idx < grid_.gdim(0);idx++)xrand[idx] = ran0() - 0.5;
    for(int idx = 0;idx < grid_.gdim(1);idx++)yrand[idx] = ran0() - 0.5;
    for(int idx = 0;idx < grid_.gdim(2);idx++)zrand[idx] = ran0() - 0.5;

    for(int ix = 0;ix < dim_[0];ix++) {

        int iix = (ix+gpt)*incx_;

        for(int iy = 0;iy < dim_[1];iy++) {

            int iiy=(iy+gpt)*incy_+iix;

            for(int iz = 0;iz < dim_[2];iz++) {
            
                int iiz=iiy+iz+gpt;
                uu_[iiz]= xrand[xoff+ix]
                        * yrand[yoff+iy]
                        * zrand[zoff+iz];
            } 
       
        } 

    } 

    updated_boundaries_=false;

    delete[] xrand;
    delete[] yrand;
    delete[] zrand;
#ifdef DEBUG
    if( mype_env().onpe0() )
        cout<<" init_rand done"<<endl;
#endif
    
} 

template <typename T>
void GridFunc<T>::initiateExchangeNorthSouth()
{
    if( directionNeumann_[1] )
    setBoundaryValuesNeumannY(valuesNeumann_[1]);

    if(mype_env().n_mpi_task(1)>1)
    {
        //cout<<" 2 PEs in direction y\n";
        const int     shift=ghost_pt();

        const size_t sizeb=shift*dim_[2]*dim_[0];
        if( north_ )
            mype_env().Irecv(&buf4_[0], sizeb, NORTH, &ns_mpireq_[3]);
        if( south_ )
            mype_env().Irecv(&buf3_[0], sizeb, SOUTH, &ns_mpireq_[1]);

        const int imax=(dim_[0]+shift)*incx_;
        const int imin=shift*incx_;
        const int jmax=shift*incy_;
        const int ymax = dim(1)*grid_.inc(1);
        const size_t sdimz=dim_[2]*sizeof(T);

        if(south_)
        {
            T* buf2_ptr = &buf2_[0];
            const T* uus=&uu_[shift*(incy_+1)];
            for (int j=0;j<jmax;j+=incy_)
            for (int i=imin; i<imax; i+=incx_ )
            {
                memcpy(buf2_ptr, &uus[i+j], sdimz);
                buf2_ptr+=dim_[2];
            }
            mype_env().Isend(&buf2_[0], sizeb, SOUTH, &ns_mpireq_[2]);
        }

        if(north_)
        {
            T* buf1_ptr = &buf1_[0];
            const T* uus=&uu_[shift+ymax];
            for (int j=0;j<jmax;j+=incy_)
            for (int i=imin; i<imax; i+=incx_ )
            {
                memcpy(buf1_ptr, &uus[i+j], sdimz);
                buf1_ptr+=dim_[2];
            }
            mype_env().Isend(&buf1_[0], sizeb, NORTH, &ns_mpireq_[0]);
        }
    }
}

template <typename T>
void GridFunc<T>::finishExchangeNorthSouth()
{
    finishExchangeNorthSouth_tm_.start();
    
    const int shift=ghost_pt();
    const int ymax = dim(1)*grid_.inc(1);
    const size_t sdimz=dim_[2]*sizeof(T);
    if(mype_env().n_mpi_task(1)>1)
    {
        const int imax=(dim_[0]+shift)*incx_;
        const int imin=shift*incx_;
        const int jmax=shift*incy_;
        
        if(south_)
            MPI_Wait(ns_mpireq_+2, MPI_STATUS_IGNORE);
        if(north_)
            MPI_Wait(ns_mpireq_, MPI_STATUS_IGNORE);
        if( north_ )
        {
            MPI_Wait(ns_mpireq_+3, MPI_STATUS_IGNORE);
            T* buf4_ptr = &buf4_[0];
            T* uus=&uu_[shift*(incy_+1)+ymax];
            for (int j=0;j<jmax;j+=incy_)
            for (int i=imin; i<imax; i+=incx_ )
            {
                memcpy(&uus[i+j], buf4_ptr, sdimz);
                buf4_ptr+=dim_[2];
            }
        }

        if( south_ )
        {
            MPI_Wait(ns_mpireq_+1, MPI_STATUS_IGNORE);
            T* buf3_ptr = &buf3_[0];
            T* uus=&uu_[shift];
            for (int j=0;j<jmax;j+=incy_)
            for (int i=imin; i<imax; i+=incx_ )
            {
                memcpy(&uus[i+j], buf3_ptr, sdimz);
                buf3_ptr+=dim_[2];
            }
        }
    }
    else // mype_env().n_mpi_task(1)==1
    {
        if(directionPeriodic_[1])
        {
            // only for i already initialized
            for (int j=0;j<shift*incy_;j+=incy_)
            for (int i=shift*incx_; i<(dim_[0]+shift)*incx_; i+=incx_ )
            {
                memcpy(&uu_[i+shift+j],               &uu_[i+shift+j+ymax],     sdimz);
                memcpy(&uu_[i+shift*(incy_+1)+j+ymax], &uu_[i+shift*(incy_+1)+j], sdimz);
     
            } 
        }
    }
    
    finishExchangeNorthSouth_tm_.stop();
}

template <typename T>
void GridFunc<T>::initiateExchangeUpDown()
{
    if( directionNeumann_[2] )
    setBoundaryValuesNeumannZ(valuesNeumann_[2]);

    const int     shift=ghost_pt();
    const int     dimxy=(dim(1)+2*ghost_pt())*dim(0);
    const int     zmax = dim(2);
    const int     iinit=(dim_[1]+2*shift)*shift;        
    const int     ione=1;

    if(mype_env().n_mpi_task(2)>1)
    {
        //int icount=0;
        const int sizeb=shift*dimxy;
        if(down_)
        {
            //icount++;
            mype_env().Irecv(&buf3_[0], sizeb, DOWN, &ud_mpireq_[1]);
        }
        if(up_)
        {
            //icount++;
            mype_env().Irecv(&buf4_[0], sizeb, UP, &ud_mpireq_[3]);
        }

        const T* const uus=&uu_[shift+incy_*iinit];
        
        if(up_){
            T* buf1_ptr=&buf1_[0];
            for(int j=0;j<shift;j++)
            {
                Tcopy(&dimxy, &uus[zmax-1-j], &incy_, buf1_ptr, &ione);
                    buf1_ptr+=dimxy;
            }
            //icount++;
            mype_env().Isend(&buf1_[0], sizeb, UP, &ud_mpireq_[0]);
        }
        if(down_)
        {
            T* buf2_ptr=&buf2_[0];
            for(int j=0;j<shift;j++)
            {
                Tcopy(&dimxy, &uus[j], &incy_, buf2_ptr, &ione);
                    buf2_ptr+=dimxy;
            }
            //icount++;
            mype_env().Isend(&buf2_[0], sizeb, DOWN, &ud_mpireq_[2]);
        }
    }
}

template <typename T>
void GridFunc<T>::finishExchangeUpDown()
{     
    finishExchangeUpDown_tm_.start();
    
    const int     dimxy=(dim(1)+2*ghost_pt())*dim(0);
    const int     zmax = dim(2);
    const int     shift=ghost_pt();
    const int     iinit=(dim_[1]+2*shift)*shift;        
    const int     ione=1;
    
    if(mype_env().n_mpi_task(2)>1)
    {
        //cout<<icount<<" communications by task "<<mype_env().mytask_()<<endl;

            //MPI_Waitall(4, ud_mpireq_, MPI_STATUS_IGNORE);
        if(up_)
                MPI_Wait(ud_mpireq_, MPI_STATUS_IGNORE);
        if(down_)
        {
                T* buf3_ptr=&buf3_[0];

                MPI_Wait(ud_mpireq_+2, MPI_STATUS_IGNORE);        
                MPI_Wait(ud_mpireq_+1, MPI_STATUS_IGNORE);
            
            for(int j=0;j<shift;j++){
                T* const uus=&uu_[shift-1-j];
                Tcopy(&dimxy, buf3_ptr, &ione, &uus[incy_*iinit], &incy_);
                        buf3_ptr+=dimxy;
            }
        }
        if(up_)
        {
            T* buf4_ptr=&buf4_[0];
                MPI_Wait(ud_mpireq_+3, MPI_STATUS_IGNORE);
            for(int j=0;j<shift;j++){
                T* const uus=&uu_[shift+zmax+j];
                Tcopy(&dimxy, buf4_ptr, &ione, &uus[incy_*iinit], &incy_);
                        buf4_ptr+=dimxy;
            }
        }
 
    }
    else
    {
        if( directionPeriodic_[2] ) /* mype_env().n_mpi_task(2)==1 */
        {
            for(int j=0;j<shift;j++)
            {
                Tcopy(&dimxy, &uu_[shift-j+iinit*incy_+zmax-1], &incy_, 
                              &uu_[shift-j+iinit*incy_-1], &incy_);
            }
            for(int j=0;j<shift;j++)
            {
                Tcopy(&dimxy, &uu_[shift+j+iinit*incy_], &incy_, 
                              &uu_[shift+j+iinit*incy_+zmax], &incy_);
            }

        }
    }
    
    finishExchangeUpDown_tm_.stop();
}

template <typename T>
void GridFunc<T>::initiateExchangeEastWest()
{
    if( directionNeumann_[0] )
        setBoundaryValuesNeumannX(valuesNeumann_[0]);

    const int     xmax = dim(0)*grid_.inc(0);
    const int size=ghost_pt()*incx_;

    if(mype_env().n_mpi_task(0)>1)
    {
        //cout<<" 2 PEs in direction x\n";
 
        /* Non-blocking MPI */
        if( east_ )
            mype_env().Irecv(&uu_[xmax+size], size, EAST, &ew_mpireq_[3]);
        if( west_ )
            mype_env().Irecv(&uu_[0],         size, WEST, &ew_mpireq_[2]);
        if( west_ )
            mype_env().Isend(&uu_[size],      size, WEST, &ew_mpireq_[1]);
        if( east_ )
            mype_env().Isend(&uu_[xmax],      size, EAST, &ew_mpireq_[0]);
    }
}

template <typename T>
void GridFunc<T>::finishExchangeEastWest()
{
    finishExchangeEastWest_tm_.start();
    
    const int xmax = dim(0)*grid_.inc(0);
    const int size=ghost_pt()*incx_;
    
    if(mype_env().n_mpi_task(0)>1)
    {
        if(west_)
            MPI_Wait(ew_mpireq_+1, MPI_STATUS_IGNORE);
        if(east_)
        {
            MPI_Wait(ew_mpireq_,   MPI_STATUS_IGNORE);
        
            MPI_Wait(ew_mpireq_+3, MPI_STATUS_IGNORE);
        }
        if(west_)
            MPI_Wait(ew_mpireq_+2, MPI_STATUS_IGNORE);
    }
    else
    {
        if( directionPeriodic_[0] )
        { /* mype_env().n_mpi_task(0)==1 */
            memcpy(&uu_[0],        &uu_[xmax],size*sizeof(T));
            memcpy(&uu_[size+xmax],&uu_[size],size*sizeof(T));
        }
    }
    
    finishExchangeEastWest_tm_.stop();
}

template <typename T>
void GridFunc<T>::setBoundaryValuesBeforeTrade()
{
    if( directionDirichlet_[0] 
     || directionDirichlet_[1] 
     || directionDirichlet_[2] )
        setBoundaryValues(0., directionDirichlet_);
    if( directionMultipole_[0] 
     || directionMultipole_[1] 
     || directionMultipole_[2] ){
        assert( bc_func_!=NULL );
        setBoundaryValues(*bc_func_, directionMultipole_);
    }
}

template <typename T>
void GridFunc<T>::defaultTrade_boundaries()
{
    if( !grid_.active() )return;
    if( updated_boundaries_ )return;
    if( ghost_pt()==0)return;

    trade_bc_tm_.start();
    
    //cout<<"GridFunc<T>::defaultTrade_boundaries(), dim_[0]="<<dim_[0]<<endl;

    setBoundaryValuesBeforeTrade();
    
    // Y direction

    initiateExchangeNorthSouth();

    finishExchangeNorthSouth();

    // Z direction

    initiateExchangeUpDown();
    
    finishExchangeUpDown();

    // X direction

    initiateExchangeEastWest();
    
    finishExchangeEastWest();

    updated_boundaries_=true;

    trade_bc_tm_.stop();
}

// set boundaries to value alpha. 
// Boundaries include ghosts points, and last layer 
// (to have an odd number of grid points with nonzero values)
template <typename T>
void GridFunc<T>::setBoundaryValues(const T alpha, const bool direction[3])
{
    //cout<<"GridFunc<T>::setBoundaryValues(const T, const bool direction[3])"<<endl;
    if( !grid_.active() )return;
    
    // early return if no direction to set BC
    if( !direction[0] && !direction[1] && !direction[2] )return;

    const int     shift= ghost_pt();
    
    assert(dim_[0]>=shift);
    assert(dim_[1]>=shift);
    assert(dim_[2]>=shift);
        
    int icount=0;
    
    const bool mympix0 = (mype_env().my_mpi(0)==0);
    const bool mympiy0 = (mype_env().my_mpi(1)==0);
    const bool mympiz0 = (mype_env().my_mpi(2)==0);
    
    const bool mympixl = (mype_env().my_mpi(0)==(mype_env().n_mpi_task(0)-1));
    const bool mympiyl = (mype_env().my_mpi(1)==(mype_env().n_mpi_task(1)-1));
    const bool mympizl = (mype_env().my_mpi(2)==(mype_env().n_mpi_task(2)-1));

    // set parameters to avoid setting values twice
    int     i0=0;
    int     j0=0;
    int     i1, j1;

    if(mympix0) i0=(shift+1)*incx_;

    if(mympixl) i1=(shift+dim_[0])*incx_;
    else        i1=(2*shift+dim_[0])*incx_;
    
    if(mympiy0) j0=(shift+1)*incy_;

    if(mympiyl) j1=(shift+dim_[1])*incy_;
    else        j1=(2*shift+dim_[1])*incy_;

    T* const pu=&uu_[0];
    
    // loop over "yz" layers
    for (int i=i0;i<i1;i+=incx_){

        // set BC in direction z
        if( direction[2] )
        for (int j=j0; j<j1; j+=incy_) {
 
            const int init_ij=i+j;

            // "low" z
            if(mympiz0){
                const int kmax=init_ij+shift;
                    for(int k=init_ij;k<=kmax;k++){
                    assert(k<grid_.sizeg());
                        pu[k]=alpha;
                    }
                icount+=shift+1;
            }

            // "high" z
            if(mympizl){
                T* const ppu = &pu[init_ij+dim_[2]+shift];
                    for(int k=0;k<shift;k++){
                    assert(k<grid_.sizeg());
                            ppu[k]=alpha;
                    }
                icount += shift;
            }
            }

        // set BC in direction y
            if( direction[1] ){
            // "low" y
            if(mympiy0){
                T* const ppu = &pu[i];
                const int kmax=incy_*(shift+1);
                for(int k=0;k<kmax;k++){
                    assert(k<grid_.sizeg());
                        ppu[k]=alpha;
                    }
                icount+=kmax;
            }


            // "high" y
            if(mympiyl){
                const int jmin=(dim_[1]+shift)*incy_;
                const int jmax=(dim_[1]+2*shift)*incy_;
                    for(int j=jmin; j<jmax; j+=incy_){
                    T* const ppu = &pu[i+j];
                    for(int k=0;k<incy_;k++){
                        assert(k<grid_.sizeg());
                        ppu[k]=alpha;
                    }
                    icount+=incy_;
                }
            }
        }

    }


    // set BC in direction x
    if( direction[0] ){
        const int size=(shift+1)*incx_;
        // 1st layers along x
        if(mympix0)
        {
            for(int i=0;i<size;i++){
                pu[i]=alpha;
            }
            icount+=size;
        }

        // last layers along x
        if(mympixl)
        {
            const int istart=incx_*(shift+dim_[0]);
            const int imax  =incx_*(2*shift+dim_[0]);
            for(int i=istart;i<imax;i++){
                assert( i<grid_.sizeg() );
                pu[i]=alpha;
            }
            icount+=(size-incx_);
        }
    }

#ifdef DEBUG
    icount-=grid_.sizeg();
    icount+=grid_.size();
    icount-=(grid_.dim(0)*grid_.dim(1));
    icount-=(grid_.dim(0)*(grid_.dim(2)-1));
    icount-=((grid_.dim(1)-1)*(grid_.dim(2)-1));

    //cout<<" icount="<<icount<<endl;
    
    if(mype_env().my_mpi(0)!=0)icount+=shift*(grid_.dim(1)-1)*(grid_.dim(2)-1);
    if(mype_env().my_mpi(1)!=0)icount+=shift*(grid_.dim(0)-1)*(grid_.dim(2)-1);
    if(mype_env().my_mpi(2)!=0)icount+=shift*(grid_.dim(1)-1)*(grid_.dim(0)-1);

    if(mype_env().my_mpi(0)!=(mype_env().n_mpi_task(0)-1))
        icount+=(shift+1)*(grid_.dim(1)-1)*(grid_.dim(2)-1);
    if(mype_env().my_mpi(1)!=(mype_env().n_mpi_task(1)-1))
        icount+=(shift+1)*(grid_.dim(0)-1)*(grid_.dim(2)-1);
    if(mype_env().my_mpi(2)!=(mype_env().n_mpi_task(2)-1))
        icount+=(shift+1)*(grid_.dim(1)-1)*(grid_.dim(0)-1);

    if((mype_env().my_mpi(0)!=0)&&(mype_env().my_mpi(1)!=0))
        icount+=shift*shift*grid_.dim(2);
    if((mype_env().my_mpi(2)!=0)&&(mype_env().my_mpi(1)!=0))
        icount+=shift*shift*grid_.dim(0);
    if((mype_env().my_mpi(0)!=0)&&(mype_env().my_mpi(2)!=0))
        icount+=shift*shift*grid_.dim(1);
    
    if(  (mype_env().my_mpi(0)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(1)!=(mype_env().n_mpi_task(1)-1)))
        icount+=(shift+1)*(shift+1)*grid_.dim(2);
    if(  (mype_env().my_mpi(0)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(2)!=(mype_env().n_mpi_task(2)-1)))
        icount+=(shift+1)*(shift+1)*grid_.dim(1);
    if(  (mype_env().my_mpi(2)!=(mype_env().n_mpi_task(2)-1))
       &&(mype_env().my_mpi(1)!=(mype_env().n_mpi_task(1)-1)))
        icount+=(shift+1)*(shift+1)*grid_.dim(0);

    if(  (mype_env().my_mpi(0)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(1)!=0))
        icount+=(shift+1)*(shift)*grid_.dim(2);
    if(  (mype_env().my_mpi(2)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(1)!=0))
        icount+=(shift+1)*(shift)*grid_.dim(0);
    if(  (mype_env().my_mpi(0)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(2)!=0))
        icount+=(shift+1)*(shift)*grid_.dim(1);


    
    if((mype_env().my_mpi(0)!=0)&&(mype_env().my_mpi(1)!=0)&&(mype_env().my_mpi(2)!=0))
        icount+=shift*shift*shift;

    if(  (mype_env().my_mpi(0)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(1)!=(mype_env().n_mpi_task(1)-1))
       &&(mype_env().my_mpi(2)!=(mype_env().n_mpi_task(2)-1)))
        icount+=(shift+1)*(shift+1)*(shift+1);

    if(  (mype_env().my_mpi(0)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(1)!=(mype_env().n_mpi_task(1)-1))
       &&(mype_env().my_mpi(2)!=0))
        icount+=(shift+1)*(shift+1)*(shift);
    if(  (mype_env().my_mpi(0)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(2)!=(mype_env().n_mpi_task(1)-1))
       &&(mype_env().my_mpi(1)!=0))
        icount+=(shift+1)*(shift+1)*(shift);
    if(  (mype_env().my_mpi(2)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(1)!=(mype_env().n_mpi_task(1)-1))
       &&(mype_env().my_mpi(0)!=0))
        icount+=(shift+1)*(shift+1)*(shift);

    if(  (mype_env().my_mpi(0)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(1)!=0)
       &&(mype_env().my_mpi(2)!=0))
        icount+=(shift+1)*(shift)*(shift);
    if(  (mype_env().my_mpi(1)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(0)!=0)
       &&(mype_env().my_mpi(2)!=0))
        icount+=(shift+1)*(shift)*(shift);
    if(  (mype_env().my_mpi(2)!=(mype_env().n_mpi_task(0)-1))
       &&(mype_env().my_mpi(1)!=0)
       &&(mype_env().my_mpi(0)!=0))
        icount+=(shift+1)*(shift)*(shift);


    //cout<<"icount="<<icount<<endl;
    //assert(icount==0);
#endif
 
}

// set boundaries to value alpha. 
// Boundaries include ghosts points, and last layer 
// (to have an odd number of grid points with nonzero values)
template <typename T>
void GridFunc<T>::setBoundaryValuesNeumannZ(const T alpha)
{
    const int     shift= ghost_pt();
    const double  hgrid=grid_.hgrid(2);

    const bool mympiz0 = (mype_env().my_mpi(2)==0);
    const bool mympizl = (mype_env().my_mpi(2)==(mype_env().n_mpi_task(2)-1));

    int     i0=0;
    int     i1=(2*shift+dim_[0])*incx_;
    int     j0=0;
    int     j1=(2*shift+dim_[1])*incy_;

    const double delta=alpha*hgrid;
    
    if( mympiz0 || mympizl )
    for (int i=i0;i<i1;i+=incx_){

        for (int j=j0; j<j1; j+=incy_) {
 
                const int init_ij=i+j;

            // "low" z
            if(mympiz0){
                const int kmax=init_ij+shift;
                    for(int k=init_ij;k<kmax;k++){
                    assert(k<grid_.sizeg());
                        uu_[k]=uu_[kmax]- (T)((double)(kmax-k)*delta);
                    }
            }

            // "high" z
            if(mympizl){
                T* const ppu = &uu_[init_ij+shift+dim_[2]-1];
                    for(int k=0;k<shift;k++){
                    assert(k<grid_.sizeg());
                            ppu[k+1]=ppu[0]+(T)((double)(k+1)*delta);
                    }
            }
            }
    }
}

template <typename T>
void GridFunc<T>::setBoundaryValuesNeumannY(const T alpha)
{
    const int     shift= ghost_pt();
    const double  hgrid=grid_.hgrid(1);
    
    const bool mympiy0 = (mype_env().my_mpi(1)==0);
    const bool mympiyl = (mype_env().my_mpi(1)==(mype_env().n_mpi_task(1)-1));

    // set parameters to avoid setting values twice
    int i0=shift*incx_;
    int i1=(shift+dim_[0])*incx_;
    
    const double delta=alpha*hgrid;
    
    if( mympiy0 || mympiyl )
    for (int i=i0;i<i1;i+=incx_){

        // set BC in direction y
        // "low" y
        if(mympiy0){
                for(int j=0; j<shift; j++){
                T* const ppu = &uu_[i+j*incy_];
                const T corrj=(T)((double)(shift-j)*delta);
                for(int k=0;k<incy_;k++){
                    assert(k<grid_.sizeg());
                        ppu[k]=uu_[i+shift*incy_+k]-corrj;
                    }
            }
        }

        // "high" y
        if(mympiyl){
            const int jmin=dim_[1]+shift;
            const int jmax=dim_[1]+2*shift;
                for(int j=jmin; j<jmax; j++){
                T* const ppu = &uu_[i+j*incy_];
                const T corrj=(T)((double)(j-jmin+1)*delta);
                for(int k=0;k<incy_;k++){
                    assert(k<grid_.sizeg());
                        ppu[k]=uu_[i+jmin*incy_-incy_+k]+corrj;
                    }
            }
        }
    }
}

template <typename T>
void GridFunc<T>::setBoundaryValuesNeumannX(const T alpha)
{
    const int     shift= ghost_pt();
    const double  hgrid=grid_.hgrid(0);
    
    const bool mympix0 = (mype_env().my_mpi(0)==0);    
    const bool mympixl = (mype_env().my_mpi(0)==(mype_env().n_mpi_task(0)-1));

    const double delta=alpha*hgrid;
    
    // 1st layers along x
    if(mympix0)
    {
        for(int i=0;i<shift;i++){
            const T corri=(T)((double)(shift-i)*delta);
            for(int k=0;k<incx_;k++)
                uu_[i*incx_+k]=uu_[shift*incx_+k]-corri;
        }
    }

    // last layers along x
    if(mympixl)
    {
        for(int i=shift+dim_[0];i<2*shift+dim_[0];i++){
            const T corri=(T)((double)(i-shift-dim_[0]+1)*delta);
            for(int k=0;k<incx_;k++)
                uu_[i*incx_+k]=uu_[(shift+dim_[0]-1)*incx_+k]+corri;
        }
    }
    
}

template <typename T>
void GridFunc<T>::setBoundaryValuesNeumann(const T alpha[3], const bool direction[3])
{
    //cout<<"GridFunc<T>::setBoundaryValuesNeumann()"<<endl;
    if( !grid_.active() )return;
    
    if( direction[1] )
        setBoundaryValuesNeumannY(alpha[1]);

    if( direction[2] )
        setBoundaryValuesNeumannZ(alpha[2]);

    if( direction[0] )
        setBoundaryValuesNeumannX(alpha[0]);
}

// set boundaries to value in GridFunc<T> values. 
// Boundaries include ghosts points, and first layer 
// (to have an odd number of grid points with nonzero values)
template <typename T>
void GridFunc<T>::setBoundaryValues(const GridFunc<T>& values, const bool direction[3])
{
    //cout<<"GridFunc<T>::setBoundaryValues(const GridFunc<T>&, const bool direction[3])"<<endl;
    assert( dim_[0]==values.dim(0) );
    assert( dim_[1]==values.dim(1) );
    assert( dim_[2]==values.dim(2) );
    assert( ghost_pt()==values.grid().ghost_pt() );
    
    if( !grid_.active() )return;
    
    // early return if no direction to set BC
    if( !direction[0] && !direction[1] && !direction[2] )return;

    const int     shift= ghost_pt();
    
    assert(dim_[0]>=shift);
    assert(dim_[1]>=shift);
    assert(dim_[2]>=shift);
        
    int icount=0;
    
    const bool mympix0 = (mype_env().my_mpi(0)==0);
    const bool mympiy0 = (mype_env().my_mpi(1)==0);
    const bool mympiz0 = (mype_env().my_mpi(2)==0);
    
    const bool mympixl = (mype_env().my_mpi(0)==(mype_env().n_mpi_task(0)-1));
    const bool mympiyl = (mype_env().my_mpi(1)==(mype_env().n_mpi_task(1)-1));
    const bool mympizl = (mype_env().my_mpi(2)==(mype_env().n_mpi_task(2)-1));

    // set parameters to avoid setting values twice
    int     i0=0;
    int     j0=0;
    int     i1, j1;

    if(mympix0) i0=(shift+1)*incx_;

    if(mympixl) i1=(shift+dim_[0])*incx_;
    else        i1=(2*shift+dim_[0])*incx_;
    
    if(mympiy0) j0=(shift+1)*incy_;

    if(mympiyl) j1=(shift+dim_[1])*incy_;
    else        j1=(2*shift+dim_[1])*incy_;

    T* const       pu=&uu_[0];
    const T* const pv=&values.uu_[0];
    
    // loop over "yz" layers
    for (int i=i0;i<i1;i+=incx_){

        // set BC in direction z
            if( direction[2] )
        for (int j=j0; j<j1; j+=incy_) {
 
                const int init_ij=i+j;

            // "low" z
            if(mympiz0){
                const int kmax=init_ij+shift;
                    for(int k=init_ij;k<=kmax;k++){
                    assert(k<grid_.sizeg());
                        pu[k]=pv[k];
                    }
                icount+=shift+1;
            }

            // "high" z
            if(mympizl){
                T* const       ppu = &pu[init_ij+dim_[2]+shift];
                const T* const ppv = &pv[init_ij+dim_[2]+shift];
                    for(int k=0;k<shift;k++){
                    assert(k<grid_.sizeg());
                            ppu[k]=ppv[k];
                    }
                icount += shift;
            }
            }

        // set BC in direction y
            if( direction[1] ){
            // "low" y
            if(mympiy0){
                T* const       ppu = &pu[i];
                const T* const ppv = &pv[i];
                const int kmax=incy_*(shift+1);
                for(int k=0;k<kmax;k++){
                    assert(k<grid_.sizeg());
                        ppu[k]=ppv[k];
                    }
                icount+=kmax;
            }


            // "high" y
            if(mympiyl){
                const int jmin=(dim_[1]+shift)*incy_;
                const int jmax=(dim_[1]+2*shift)*incy_;
                    for(int j=jmin; j<jmax; j+=incy_){
                    T* const       ppu = &pu[i+j];
                    const T* const ppv = &pv[i+j];
                    for(int k=0;k<incy_;k++){
                        assert(k<grid_.sizeg());
                            ppu[k]=ppv[k];
                        }
                    icount+=incy_;
                }
            }
        }

    }


    // set BC in direction x
    if( direction[0] ){
        const int size=(shift+1)*incx_;
        // 1st layers along x
        if(mympix0)
        {
            for(int i=0;i<size;i++){
                pu[i]=pv[i];
            }
            icount+=size;
        }

        // last layers along x
        if(mympixl)
        {
            const int istart=incx_*(shift+dim_[0]);
            const int imax  =incx_*(2*shift+dim_[0]);
            for(int i=istart;i<imax;i++){
                assert( i<grid_.sizeg() );
                pu[i]=pv[i];
            }
            icount+=(size-incx_);
        }
    }

}

// dot product on the global grid (distributed)
/* This is split into the double-type argument and float-type argument 
 * below. This is necessary to ensure that the underlying MPdot routine 
 * that is called returns the same result whether T=float and vv is double
 * or T=double and vv is float. Otherwise the GridFunc copy constructor 
 * would use a copy of vv of the same type as "this" GridFunc object (i.e. type T).
 * This could lead to different results for float-double and double-float combinations 
 * of this function.
*/
template <typename T>
double GridFunc<T>::gdot(const GridFunc<double>& vv)const
{
    if( !grid_.active() )return 0.;

    const int     nghosts=ghost_pt();

    int           dimz=grid_.dim(2);
    
    const int endx =(nghosts+dim_[0])*incx_;
    const int endy =(nghosts+dim_[1])*incy_;

    int initx = nghosts*incx_;
    int inity = nghosts*incy_;
    int initz = nghosts;
   
    // remove "layers" belonging to BC
    if( ((bc_[0]!=1)||(vv.bc(0)!=1)) && mype_env().my_mpi(0)==0 ) initx+=incx_;
    if( ((bc_[1]!=1)||(vv.bc(1)!=1)) && mype_env().my_mpi(1)==0 ) inity+=incy_;
    if( ((bc_[2]!=1)||(vv.bc(2)!=1)) && mype_env().my_mpi(2)==0 ){ 
        dimz--;
        initz++;
    }
    
    const double * const vv1=vv.uu();
    const T * const vv2=uu_;
    double  my_dot=0.;
    for(int ix=initx;ix<endx;ix+=incx_){
        for(int iy=inity;iy<endy;iy+=incy_){
            int iz=ix+iy+initz;
            my_dot+=MPdot(dimz, &vv1[iz], &vv2[iz]);
        }
    }

#ifdef USE_MPI
    MGmol_MPI& mmpi=*(MGmol_MPI::instance());
    if( mype_env().n_mpi_tasks()>1 ){
        double  sum=0.;
        int rc=mmpi.allreduce(&my_dot, &sum, 1, MPI_SUM);
        if(rc!=MPI_SUCCESS){
            cout<<"MPI_Allreduce double sum failed in gdot!!!"<<endl;
            mype_env().globalExit(2);
        }
        my_dot=sum;
    }
#endif

    return my_dot;
    
}

// dot product on the global grid (distributed)
template <typename T>
double GridFunc<T>::gdot(const GridFunc<float>& vv)const
{
    if( !grid_.active() )return 0.;

    const int     nghosts=ghost_pt();

    int           dimz=grid_.dim(2);
    
    const int endx =(nghosts+dim_[0])*incx_;
    const int endy =(nghosts+dim_[1])*incy_;

    int initx = nghosts*incx_;
    int inity = nghosts*incy_;
    int initz = nghosts;
   
    // remove "layers" belonging to BC
    if( ((bc_[0]!=1)||(vv.bc(0)!=1)) && mype_env().my_mpi(0)==0 ) initx+=incx_;
    if( ((bc_[1]!=1)||(vv.bc(1)!=1)) && mype_env().my_mpi(1)==0 ) inity+=incy_;
    if( ((bc_[2]!=1)||(vv.bc(2)!=1)) && mype_env().my_mpi(2)==0 ){ 
        dimz--;
        initz++;
    }
    
    const float * const vv1=vv.uu();
    const T * const vv2=uu_;
    double  my_dot=0.;
    for(int ix=initx;ix<endx;ix+=incx_){
        for(int iy=inity;iy<endy;iy+=incy_){
            int iz=ix+iy+initz;
            my_dot+=MPdot(dimz, &vv1[iz], &vv2[iz]);
        }
    }

#ifdef USE_MPI
    MGmol_MPI& mmpi=*(MGmol_MPI::instance());
    if( mype_env().n_mpi_tasks()>1 ){
        double  sum=0.;
        int rc=mmpi.allreduce(&my_dot, &sum, 1, MPI_SUM);
        if(rc!=MPI_SUCCESS){
            cout<<"MPI_Allreduce double sum failed in gdot!!!"<<endl;
            mype_env().globalExit(2);
        }
        my_dot=sum;
    }
#endif

    return my_dot;
    
}

template <typename T>
double GridFunc<T>::norm2()const
{
    double my_dot=gdot(*this);

    return sqrt(my_dot*grid_.vel());
    
}

template <typename T>
bool GridFunc<T>::def_const()const
{
    if( !grid_.active() )return false;

    int  tmp=0;
    int  sum=0;
    if(( !bc_[0] || !bc_[1] || !bc_[2] ))tmp=1;

    MGmol_MPI& mmpi=*(MGmol_MPI::instance());
    if( mype_env().n_mpi_tasks()>1 )
    {
        int rc=mmpi.allreduce(&tmp, &sum, 1, MPI_SUM);
        if(rc!=MPI_SUCCESS){
            cout<<"MPI_Allreduce double sum failed!!!"<<endl;
            mype_env().globalExit(2);
        }
    }

    if(sum>0)return true;
    else     return false;
}

template <typename T>
double GridFunc<T>::get_average()
{
    if( !grid_.active() )return 0.;

    double  sum=0.;
    const int gpt=ghost_pt();
 
    assert( bc_[0]==1 && bc_[1]==1 && bc_[2]==1 );

    for(int ix = gpt;ix < gpt+dim_[0];ix++)
    {
        int iix = ix*incx_;

        for(int iy = gpt;iy < gpt+dim_[1];iy++)
        {
            int iiy=iy*incy_+iix;
            const T* __restrict__ pu=&uu_[iiy+gpt];
            
            for(int iz = 0;iz < dim_[2];iz++)
            {            
                sum+=(double)pu[iz];
            }
        } 

    }    
    sum/=(double)grid_.size();
    //cout<<"sum="<<sum<<endl;

    MGmol_MPI& mmpi=*(MGmol_MPI::instance());
    if( mype_env().n_mpi_tasks()>1 )
    {
        double  tmp=0.;
        int rc=mmpi.allreduce(&sum, &tmp, 1, MPI_SUM);
        if(rc!=MPI_SUCCESS){
            cout<<"MPI_Allreduce double sum failed!!!"<<endl;
            mype_env().globalExit(2);
        }
        sum=tmp/((double)mype_env().n_mpi_tasks());
    }
    
    return sum;
}

template <typename T>
double GridFunc<T>::average0()
{
    const double sum=get_average();

    //cout<<"average to 0.\n";
    if(!def_const()){
#ifdef DEBUG
        cout<<"substract "<<sum<<" to function"<<endl;
#endif
        *this-=sum;
    }
    return sum;
}

//check if function is zero everywhere
//(to be used in tests)
template <typename T>
bool GridFunc<T>::isZero(const double tol, const bool wghosts)
{
    if( !grid_.active() )return true;
    
    const short nghosts = ghost_pt();

    int init[3];
    int end[3];
    if( wghosts ){
        for(int i=0;i<3;i++){
            init[i]=0;
            end[i] =grid_.dim(i)+2*nghosts;
        }
    }
    else
        for(int i=0;i<3;i++){
            init[i]=nghosts;
            end[i] =grid_.dim(i)+nghosts;
        }
    
    for(int ix = init[0];ix < end[0];ix++) {

        const int iix = ix*grid_.inc(0);

        for(int iy = init[1];iy < end[1];iy++) {

            const int iiy=iy*grid_.inc(1)+iix;

            for(int iz = init[2];iz < end[2];iz++) {
            
                const int iiz=iiy+iz;

                if( fabs(uu_[iiz])>tol ){
                    cout<<"i="<<ix<<",j="<<iy<<",k="<<iz<<", u="<<uu_[iiz]<<endl;
                    return false;
                }
            } 
       
        }

    } 
    return true;
}

template <typename T>
double GridFunc<T>::integral()const
{
    if( !grid_.active() )return 0.;

    double integral=0.;
    const int gpt=ghost_pt();

    for(int ix = gpt;ix < gpt+dim_[0];ix++) {

        int iix = ix*incx_;

        for(int iy = gpt;iy < gpt+dim_[1];iy++) {

            int iiy=iy*incy_+iix;
            const T* __restrict__ pu=&uu_[iiy+gpt];
            
            for(int iz = 0;iz < dim_[2];iz++) {            
                integral+=(double)pu[iz];
            } 
       
        } 

    }    
    //cout<<"integral="<<integral<<endl;

    MGmol_MPI& mmpi=*(MGmol_MPI::instance());
    if(grid_. mype_env().n_mpi_tasks()>1 ){
        double  tmp=0.;
        int rc=mmpi.allreduce(&integral, &tmp, 1, MPI_SUM);
        if(rc!=MPI_SUCCESS){
            cout<<"MPI_Allreduce double sum failed!!!"<<endl;
            mype_env().globalExit(2);
        }
        integral=tmp;
    }
    
    return integral*grid_.vel();
}

template <typename T>
void GridFunc<T>::extend3D(GridFunc<T>& ucoarse)
{
    extend3D_tm_.start();
    
    if( !grid_.active() )return;

    //assert(bc_==ucoarse.bc());
    assert(mype_env().n_mpi_tasks()==ucoarse.grid().mype_env().n_mpi_tasks());
    assert( dim(0)>=2*ucoarse.dim(0) );
    assert( dim(1)>=2*ucoarse.dim(1) );
    assert( dim(2)>=2*ucoarse.dim(2) );
    assert( ghost_pt()>0 );

    const int incx_fine   = grid_.inc(0);
    const int incy_fine   = grid_.inc(1);
    const int incy_coarse = ucoarse.grid_.inc(1);
    const int incx_coarse = ucoarse.grid_.inc(0);

    const short nghosts_fine   = ghost_pt();
    const short nghosts_coarse = ucoarse.ghost_pt();

    const int cdimx=ucoarse.dim(0);
    const int cdimy=ucoarse.dim(1);
    const int cdimz=ucoarse.dim(2);
    
    int ione=1;
    int itwo=2;
    
    bc_[0]=ucoarse.bc_[0];
    bc_[1]=ucoarse.bc_[1];
    bc_[2]=ucoarse.bc_[2];
    

    if(!ucoarse.updated_boundaries()) ucoarse.trade_boundaries();

    // loop over coarse grid
    // to inject coarse values on fine grid, including first ghost
    for(int ix=0;ix<=cdimx;ix++){
    
        const int ifx = incx_fine*(2*ix+nghosts_fine)+nghosts_fine;
        const int icx = incx_coarse*(ix+nghosts_coarse)+nghosts_coarse;

        for(int iy=0;iy<=cdimy;iy++){

            const int ify = ifx+incy_fine*(2*iy+nghosts_fine);
            const int icy = icx+incy_coarse*(iy+nghosts_coarse);

            int size=cdimz+1;
            Tcopy(&size, &ucoarse.uu_[icy], &ione, &uu_[ify], &itwo);
        }

    }
    
    // now work with fine level only
    for(int ix=nghosts_fine;ix<dim_[0]+nghosts_fine;ix=ix+2){
    
        int ifx = incx_fine*ix;

        for(int iy=nghosts_fine;iy<dim_[1]+nghosts_fine;iy=iy+2){
        
            int ify = ifx+incy_fine*iy;

            for(int iz=nghosts_fine+1;iz<dim_[2]+nghosts_fine;iz=iz+2){
            
                int izf = ify+iz;

                uu_[izf] = 0.5 * (uu_[izf-1] + uu_[izf+1]);
                assert(izf<grid_.sizeg());
            }

        }

        for(int iy=nghosts_fine+1;iy<dim_[1]+nghosts_fine;iy=iy+2){
        
            int ify = ifx+incy_fine*iy;

            for(int iz=nghosts_fine;iz<dim_[2]+nghosts_fine;iz=iz+2){ 
            
                int izf = ify+iz;

                    uu_[izf] = 0.5 * (uu_[izf+incy_fine] + uu_[izf-incy_fine]);

                assert(izf<grid_.sizeg());
            }

            for(int iz=nghosts_fine+1;iz<dim_[2]+nghosts_fine;iz=iz+2){ 
            
                int izf = ify+iz;

                    uu_[izf] = 0.25 * (uu_[izf+1+incy_fine] +
                                            uu_[izf+1-incy_fine] +
                                            uu_[izf-1+incy_fine] +
                                            uu_[izf-1-incy_fine]);

                assert(izf<grid_.sizeg());
            }

        }

    }


    for(int ix=nghosts_fine+1;ix<dim_[0]+nghosts_fine;ix=ix+2){
    
        int ifx = incx_fine*ix;

        for(int iy=nghosts_fine;iy<dim_[1]+nghosts_fine;iy=iy+2){
        
            int ify = ifx+incy_fine*iy;

            for(int iz=nghosts_fine+1;iz<dim_[2]+nghosts_fine;iz=iz+2){ 
            
                int izf = ify+iz;

                    uu_[izf] = 0.25 *(uu_[izf+incx_fine+1] +
                                           uu_[izf+incx_fine-1] +
                                           uu_[izf-incx_fine+1] +
                                           uu_[izf-incx_fine-1]);

            }

            for(int iz=nghosts_fine;iz<dim_[2]+nghosts_fine;iz=iz+2){ 
            
                int izf = ify+iz;

                    uu_[izf] = 0.5 * (uu_[izf+incx_fine] + uu_[izf-incx_fine]);

            }

        }

        for(int iy=nghosts_fine+1;iy<dim_[1]+nghosts_fine;iy=iy+2){
        
            int ify = ifx+incy_fine*iy;

            for(int iz=nghosts_fine+1;iz<dim_[2]+nghosts_fine;iz=iz+2){ 
            
                int izf = ify+iz;

                    uu_[izf] = 0.125 *(uu_[izf+incx_fine+incy_fine+1] +
                                            uu_[izf+incx_fine+incy_fine-1] +
                                            uu_[izf+incx_fine-incy_fine+1] +
                                            uu_[izf+incx_fine-incy_fine-1] +
                                            uu_[izf-incx_fine+incy_fine+1] +
                                            uu_[izf-incx_fine+incy_fine-1] +
                                            uu_[izf-incx_fine-incy_fine+1] +
                                            uu_[izf-incx_fine-incy_fine-1]);

            }

            for(int iz=nghosts_fine;iz<dim_[2]+nghosts_fine;iz=iz+2){ 
            
                int izf = ify+iz;

                    uu_[izf] = 0.25 *(uu_[izf+incx_fine+incy_fine] +
                                           uu_[izf+incx_fine-incy_fine] +
                                           uu_[izf-incx_fine+incy_fine] +
                                           uu_[izf-incx_fine-incy_fine]);

            }

        }

    }

    updated_boundaries_=false;

    extend3D_tm_.stop();
}

template <typename T>
void GridFunc<T>::restrict3D(GridFunc<T>& ucoarse)
{
    //cout<<"GridFunc<T>::restrict3D()"<<endl;
    restrict3D_tm_.start();

    if( !grid_.active() )return;

    assert( dim(0)==2*ucoarse.dim(0) );
    assert( dim(1)==2*ucoarse.dim(1) );
    assert( dim(2)==2*ucoarse.dim(2) );

    ucoarse.set_bc(bc_[0],bc_[1],bc_[2]);

    const int incx_fine=grid_.inc(0);
    const int incy_fine=grid_.inc(1);
    const int incy_coarse = ucoarse.inc(1);
    const int incx_coarse = ucoarse.inc(0);

    const short nghosts_fine   = ghost_pt();
    const short nghosts_coarse = ucoarse.ghost_pt();
    
    if(!updated_boundaries_) trade_boundaries();
    
    int     ixf = nghosts_fine  *incx_fine;
    int     ixc = nghosts_coarse*incx_coarse;
    
    const int cdim0=ucoarse.dim(0);
    const int cdim1=ucoarse.dim(1);
    const int cdim2=ucoarse.dim(2);

    // loop over coarse grid
    for(int ix=nghosts_fine;ix<cdim0+nghosts_fine;ix++){

        int iyf = ixf+incy_fine*nghosts_fine;
        int iyc = ixc+incy_coarse*nghosts_coarse;

        for(int iy=nghosts_fine;iy<cdim1+nghosts_fine;iy++){
        
            int izf = iyf+nghosts_fine;

            const T* const u0   =uu_+izf;
            const T* const umx  =u0-incx_fine;
            const T* const upx  =u0+incx_fine;
            const T* const umy  =u0-incy_fine;
            const T* const upy  =u0+incy_fine;
            const T* const umxpy  =u0-incx_fine+incy_fine;
            const T* const upxpy  =u0+incx_fine+incy_fine;
            const T* const umxmy  =u0-incx_fine-incy_fine;
            const T* const upxmy  =u0+incx_fine-incy_fine;
            for(int iz=0;iz<cdim2;iz++){
                            
                const int twoiz=2*iz;
                
                double face =
                         (double)upx[twoiz]
                       + (double)umx[twoiz]
                       + (double)upy[twoiz]
                       + (double)umy[twoiz]
                       + (double)u0[twoiz-1]
                       + (double)u0[twoiz+1];
                
                double corner = 
                         (double)upxpy[twoiz-1] +
                         (double)upxpy[twoiz+1] +
                         (double)upxmy[twoiz-1] +
                         (double)upxmy[twoiz+1] +
                         (double)umxpy[twoiz-1] +
                         (double)umxpy[twoiz+1] +
                         (double)umxmy[twoiz-1] +
                         (double)umxmy[twoiz+1];
                                                    
                double edge   =
                         (double)upy[twoiz-1]
                       + (double)upy[twoiz+1]
                       + (double)umy[twoiz-1]
                       + (double)umy[twoiz+1]
                       + (double)upx[twoiz-1]
                       + (double)upx[twoiz+1]
                       + (double)umx[twoiz-1]
                       + (double)umx[twoiz+1]
                       + (double)umxmy[twoiz]
                       + (double)upxmy[twoiz]
                       + (double)umxpy[twoiz]
                       + (double)upxpy[twoiz];
                

                ucoarse.uu_[iyc+iz+nghosts_coarse] 
                    = (T)(inv64 * ( 8. * u0[twoiz] 
                              + 4. * face + 2. * edge + corner));

            }
            
            iyf += 2*incy_fine;
            iyc += incy_coarse;
        }
        
        ixf += 2*incx_fine;
        ixc += incx_coarse;
    }

    ucoarse.set_updated_boundaries(false);

    restrict3D_tm_.stop();
}

template <typename T>
void GridFunc<T>::test_setBoundaryValues()
{
    cout<<" test_setBoundariesValues() on grid "<<dim(0)<<","<<dim(1)<<","<<dim(2)<<endl;

    *this=1.;
 
    double center[3];
    for(int ii=0;ii<3;ii++)
        center[ii]=(T)(0.5*(grid_.ll(ii)-grid_.hgrid(ii)));

    init_radial(center,0.); 

    bool direction[3]={true,true,true};
    setBoundaryValues(0., direction);   
    
    double my_dot=gdot(*this);

    cout<<" gdim="<<grid_.gdim(0)<<","<<grid_.gdim(1)<<","<<grid_.gdim(2)<<endl;

    cout<<" dot="<<my_dot<<endl;
    cout<<" grid="<<((grid_.gsize()
           -grid_.gdim(0)*grid_.gdim(1)
           -grid_.gdim(0)*(grid_.gdim(2)-1)
           -(grid_.gdim(1)-1)*(grid_.gdim(2)-1)))<<endl;
        
    assert( fabs(my_dot-
            (grid_.gsize()
           -grid_.gdim(0)*grid_.gdim(1)
           -grid_.gdim(0)*(grid_.gdim(2)-1)
           -(grid_.gdim(1)-1)*(grid_.gdim(2)-1)) )<1.e-8 );
    

    if( (2*(dim(0)>>1)==dim(0))&&(2*(dim(1)>>1)==dim(1))&&(2*(dim(2)>>1)==dim(2)) )
    if((dim(0)>1)&&(dim(1)>1)&&(dim(2)>1))
    {
        Grid coarse_G=grid_.coarse_grid();
        GridFunc<T>  ucoarse(coarse_G,bc_[0],bc_[1],bc_[2]);

        if((coarse_G.dim(0)>=ghost_pt())&&(coarse_G.dim(1)>=ghost_pt())
           &&(coarse_G.dim(2)>=ghost_pt()))
                ucoarse.test_setBoundaryValues(); 
    }

}

template <typename T>
void GridFunc<T>::test_grid_transfer()
{
    Grid coarse_G=grid_.coarse_grid();
 
    GridFunc<T>  rcoarse(coarse_G,bc_[0],bc_[1],bc_[2]);

    if(coarse_G.dim(0)<ghost_pt())return;
    if(coarse_G.dim(1)<ghost_pt())return;
    if(coarse_G.dim(2)<ghost_pt())return;

    if( mype_env().onpe0() )
        cout<<" Test grid transfer() between grids "<<dim(0)
            <<" and "<<coarse_G.dim(0)<<endl;

    if( bc_[0]!=1 || bc_[1]!=1 || bc_[2]!=1){
        if( mype_env().onpe0() )
            cout<<" test_grid_transfer() requires periodic BC. skipped\n";
        return;
    }
    *this=2.11;
    restrict3D(rcoarse);
    rcoarse.trade_boundaries();
    memset(uu_,0,grid_.sizeg()*sizeof(T));
    extend3D(rcoarse);
    
    for(int ix = ghost_pt();ix < ghost_pt()+dim(0);ix++) {

        int iix = ix*grid_.inc(0);

        for(int iy = ghost_pt();iy < ghost_pt()+dim(1);iy++) {

            int iiy=iy*grid_.inc(1)+iix;

            for(int iz = ghost_pt();iz < ghost_pt()+dim(2);iz++) {

                if(fabs(uu_[iiy+iz]-2.11)>0.000001){
                    cout<<"Test grid transfer: u["<<ix<<" "<<iy<<" "<<iz<<"]="<<uu_[iiy+iz]<<endl;

                    std::exit(0);
                }

            }
        }
    } 
    
    if( mype_env().onpe0() )
        cout<<" Testgrid transfer() passed"<<endl;
       
    if( (2*(coarse_G.dim(0)>>1)==coarse_G.dim(0))
      &&(2*(coarse_G.dim(1)>>1)==coarse_G.dim(1))
      &&(2*(coarse_G.dim(2)>>1)==coarse_G.dim(2)) )
    if((coarse_G.dim(0)>ghost_pt())&&(coarse_G.dim(1)>ghost_pt())
        &&(coarse_G.dim(2)>ghost_pt())){

            rcoarse.test_grid_transfer();
 
    }


}

template <typename T>
void GridFunc<T>::test_newgrid()
{
    Grid    new_grid(grid_,0);
    assert( fabs(grid_.vel()-new_grid.vel())<1.e-15 );
    assert(dim(0)==dim(0));
    GridFunc<T> f(new_grid,bc_[0],bc_[1],bc_[2]);
    f.init_rand();
    
    GridFunc<T> g(f.uu_,grid_,bc_[0],bc_[1],bc_[2],'l');
    
    double s1=norm(f);
    double s2=norm(g);
 
    cout<<" norm 1="<<s1<<", norm 2="<<s2<<endl;

    assert( fabs(s1-s2)<1.e-15 );

    
    cout<<" 1st Test OK "<<endl;

    int dims=dim(0)*dim(1)*dim(2);
    T* hu=new T[dims];
    for(int ix = 0;ix < dims;ix++) {

        hu[ix]= (T)(((double)std::rand())/32768.);
        //hu[ix]=1.;

    } 
    GridFunc<T> h(hu,grid_,bc_[0],bc_[1],bc_[2],'l');


    s1=norm(h);
    cout<<" h(hu): bc="<<h.bc(0)<<","<<h.bc(1)<<","<<h.bc(2)<<endl;
  
    double  my_dot=MPdot(dims, hu, hu);

#ifdef USE_MPI
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    double  sum=0.;
    int rc=mmpi.allreduce(&my_dot, &sum, 1, MPI_SUM);
    if(rc!=MPI_SUCCESS){
        cout<<"MPI_Allreduce double sum failed!!!"<<endl;
        mype_env().globalExit(2);
    }
    my_dot=sum;
#endif
    s2=sqrt(grid_.vel()*my_dot);
 
    //cout<<" norm 1="<<s1<<", norm 2="<<s2<<endl;
    //f.print_radial("s1.dat");
    //g.print_radial("s2.dat");

    assert(grid_.inc(2)==1);
 
    //h.print("h.dat");
 
    // print hu
    const int incx1=new_grid.inc(0);
    const int incy1=new_grid.inc(1);
    ofstream  tfile("hu.dat");
    for(int ix = 0;ix < dim(0);ix++) {

       int ix1 = ix*incx1;

       for(int iy = 0;iy < dim(1);iy++) {

           int iy1 = ix1+iy*incy1;
           
           for(int iz =0;iz < dim(2);iz++) {
    
               int iz1 = iy1+iz;

               tfile<<hu[iz1];

           } 
           
           tfile<<endl;
       } 

    } 

    assert( fabs(s1-s2)<1.e-15 );
    delete[] hu;

    cout<<" 2nd Test OK "<<endl;
}

template <typename T>
void GridFunc<T>::test_trade_boundaries()
{
    int     ix, iy, iz, iix, iiy;

    cout<<" Test trade_boundaries() on grid "
        <<dim(0)<<","<<dim(1)<<","<<dim(2)<<", bc="<<bc(0)<<","<<bc(1)<<","<<bc(2)<<endl;


#ifdef USE_MPI
    int ntasks;
    MPI_Comm_size(mype_env().comm(), &ntasks);
    assert(ntasks==mype_env().n_mpi_tasks());
#endif

    memset(uu_,0,grid_.sizeg()*sizeof(T));

    for(ix = ghost_pt();ix < ghost_pt()+dim(0);ix++) {

        iix = ix*grid_.inc(0);

        for(iy = ghost_pt();iy < ghost_pt()+dim(1);iy++) {

            iiy=iy*grid_.inc(1)+iix;

            for(iz = ghost_pt();iz < ghost_pt()+dim(2);iz++) {

                uu_[iiy+iz]=1.111;

            }
        }
    }

#ifdef DEBUG
    double norm_before=norm2();
#endif

    trade_boundaries();

#ifdef DEBUG
    double norm_after=norm2();
    assert(norm_before==norm_after);
#endif

    int     initx, inity, initz, endx, endy, endz;

    if(mype_env().my_mpi(0)==0 && !bc_[0])initx=ghost_pt();
    else                                        initx=0;
    if(mype_env().my_mpi(1)==0 && !bc_[1])inity=ghost_pt();
    else                                        inity=0;
    if(mype_env().my_mpi(2)==0 && !bc_[2])initz=ghost_pt();
    else                                        initz=0;

    if(mype_env().my_mpi(0)==(mype_env().n_mpi_task(0)-1))endx=ghost_pt()+dim(0)-1;
    else                                                  endx=2*ghost_pt()+dim(0);
    if(mype_env().my_mpi(1)==(mype_env().n_mpi_task(1)-1))endy=ghost_pt()+dim(1)-1;
    else                                                  endy=2*ghost_pt()+dim(1);
    if(mype_env().my_mpi(2)==(mype_env().n_mpi_task(2)-1))endz=ghost_pt()+dim(2)-1;
    else                                                  endz=2*ghost_pt()+dim(2);



    //cout<<" u["<<0<<" "<<0<<" "<<0<<"]="<<uu_[0]<<endl;
    //cout<<" initx="<<initx<<endl;

    for(ix = initx;ix < endx;ix++) {

        iix = ix*grid_.inc(0);

        for(iy = inity;iy < endy;iy++) {

            iiy=iy*grid_.inc(1)+iix;

            for(iz = initz;iz < endz;iz++) {

                if(fabs(uu_[iiy+iz]-1.111)>0.000001){
                    cout<<" u["<<ix<<" "<<iy<<" "<<iz<<"]="<<uu_[iiy+iz]<<endl;
                    cout<<"end="<<endx<<","<<endy<<","<<endz<<endl;
                    cout<<"init="<<initx<<","<<inity<<","<<initz<<endl;

                    std::exit(0);
                }

            }
        }
    }    

    if((dim(0)>1)&&(dim(1)>1)&&(dim(2)>1))
    if( (2*(dim(0)>>1)==dim(0))&&(2*(dim(1)>>1)==dim(1))&&(2*(dim(2)>>1)==dim(2)) ){

            Grid coarse_G=grid_.coarse_grid();
        GridFunc<T>  ucoarse(coarse_G,bc_[0],bc_[1],bc_[2]);

        if((coarse_G.dim(0)>ghost_pt())&&(coarse_G.dim(1)>ghost_pt())
           &&(coarse_G.dim(2)>ghost_pt()))
                ucoarse.test_trade_boundaries();
 
    }

}

template <typename T>
void GridFunc<T>::jacobi(const GridFunc<T>& v, const GridFunc<T>& epsilon, const double c0)
{
    if( !grid_.active() )return;

    T* __restrict__ pu=&uu_[0];
    T* __restrict__ vv=&v.uu_[0];
    T* __restrict__ ee=&epsilon.uu_[0];
    const int sizeg=grid().sizeg();

    for(int j=0;j<sizeg;j++)
    {
#ifdef DEBUG
        if(ee[j]<=0.)cout<<"ee["<<j<<"]="<<ee[j]<<endl;
        assert(ee[j]>0.);
#endif
        pu[j] -= (T)((1./(c0*ee[j]))*vv[j]);
    }

}

template <typename T>
void GridFunc<T>::add_prod(const GridFunc<T>& v1, const GridFunc<T>& v2)
{
    if( !grid_.active() )return;

    assert(v1.grid_.sizeg()==grid_.sizeg());
    assert(v2.grid_.sizeg()==grid_.sizeg());

    T* __restrict__ pu=&uu_[0];
    const T* __restrict__ vv1=&v1.uu_[0];
    const T* __restrict__ vv2=&v2.uu_[0];
    const int sizeg=grid_.sizeg();

    for(int i=0;i<sizeg;i++)
    {
        pu[i]+= (T)(vv1[i]*vv2[i]);
    }

    updated_boundaries_=false;
}

template <typename T>
void GridFunc<T>::substract_prod(const GridFunc<T>& v1, const GridFunc<T>& v2)
{
    if( !grid_.active() )return;

    assert(v1.grid_.sizeg()==grid_.sizeg());
    assert(v2.grid_.sizeg()==grid_.sizeg());

    T* __restrict__ pu=&uu_[0];
    const T* __restrict__ vv1=&v1.uu_[0];
    const T* __restrict__ vv2=&v2.uu_[0];
    const int sizeg=grid_.sizeg();

    //cout<<"begin substract...\n";
    for(int i=0;i<sizeg;i++)
    {
        pu[i]-= (T)(vv1[i]*vv2[i]);
    }

    //cout<<"end substract...\n"<<fflush;

    updated_boundaries_=false;
}

template <typename T>
void GridFunc<T>::smooth_by_coarsening(int nlevel)
{
    if( !grid_.active() )return;

    if(nlevel <= 0){
        return;
    }

    Grid      coarse_grid=grid().coarse_grid();
    GridFunc<T>  gf_coarse(coarse_grid,bc_[0],bc_[1],bc_[2]);
    restrict3D(gf_coarse);
    
    gf_coarse.smooth_by_coarsening(nlevel-1);
    
    extend3D(gf_coarse);
    
}

template <typename T>
void GridFunc<T>::add_bias(const double bias)
{
    if( !grid_.active() )return;

    if( bc_[0] && bc_[1] && bc_[2]){
        if(mype_env().onpe0())
            cout<<"Add bias "<<bias<<endl;

        const short shift=ghost_pt();
        int istart=grid_.dim(0)*mype_env().my_mpi(0);
        for(int i=0;i<grid_.dim(0);i++)
        for(int j=0;j<grid_.dim(1);j++)
        for(int k=0;k<grid_.dim(2);k++)
            uu_[incx_*(i+shift)+incy_*(j+shift)+shift+k]
                -= (T)(((double)(i+istart)/(double)(grid_.gdim(0)))*bias);    
        average0();
    }
}

template <typename T>
double GridFunc<T>::get_bias()
{
    if( !grid_.active() )return 0.;

    double  bias=0;
    if(bc_[0] && bc_[1] && bc_[2]){
        //cout<<"Compute bias "<<bias<<endl;

        const short shift=ghost_pt();
        int istart=grid_.dim(0)*mype_env().my_mpi(0);
        int i0=1;
        if(istart>0)i0=0;

        T* ref;
        if(mytask_==0){
            ref=uu()+shift*incx_;
        }else{
            ref=new T[incx_];
        }
#ifdef USE_MPI
        MGmol_MPI& mmpi=*(MGmol_MPI::instance());
        int rc=mmpi.bcast(ref, incx_, 0);
        if(rc!=MPI_SUCCESS){
            cout<<"MPI_Bcast failed in get_bias()!!!"<<endl;
            mype_env().globalExit(2);
        }
#endif

        for(int j=0;j<grid_.dim(1);j++)
        for(int k=0;k<grid_.dim(2);k++){
            for(int i=i0;i<grid_.dim(0);i++){
                assert((i+istart)>0);
                bias-= (uu_[incx_*(i+shift)+incy_*(j+shift)+shift+k]
                       -ref[incy_*(j+shift)+shift+k])
                         /((double)(i+istart)/(double)(grid_.gdim(0)));    
            }
        }
        bias=mype_env().double_sum_all(bias);
        //bias /= (grid_.gsize()
        //        -grid_.gdim(1)*grid_.gdim(2));
        bias /= (grid_.gdim(1)*grid_.gdim(2));
        bias /= (grid_.gdim(0)-1);
        if(mytask_!=0)
            delete[] ref;
    }
    return bias;
}

template <typename T>
GridFunc<T>::~GridFunc<T>()
{
    if( uu_!=0 )
        delete[] uu_;

    //cout<<"Destructor for function on grid "<<grid_.level()<<endl;    
}

template <typename T>
void GridFunc<T>::getCellCornersValues(const int i, const int j, const int k,
                                    double val[8])const
{
    const short shift=ghost_pt();
    T* pu=&uu_[shift*incx_+shift*incy_+shift];
    val[0]=(double)pu[    i*incx_+    j*incy_+k  ];
    val[1]=(double)pu[(i+1)*incx_+    j*incy_+k  ];
    val[2]=(double)pu[    i*incx_+(j+1)*incy_+k  ];
    val[3]=(double)pu[(i+1)*incx_+(j+1)*incy_+k  ];
    val[4]=(double)pu[    i*incx_+    j*incy_+k+1];
    val[5]=(double)pu[(i+1)*incx_+    j*incy_+k+1];
    val[6]=(double)pu[    i*incx_+(j+1)*incy_+k+1];
    val[7]=(double)pu[(i+1)*incx_+(j+1)*incy_+k+1];
}

// assign values from other GridFunc<T> possibly "global" (replicated)
template <typename T>
void GridFunc<T>::assign(const GridFunc<T>& src, const char dis)
{
    assert(dis=='d' || dis=='g' );

    if( !grid_.active() )return;

    if( (&grid_ == &src.grid_) && dis=='d' )
    {
        memcpy(uu_, src.uu_, grid_.sizeg()*sizeof(T));
    }
    else //more expensive assign
    {
    
    const int incx_dst   =grid_.inc(0);
    const int incy_dst   =grid_.inc(1);
    const short nghosts_dst=ghost_pt();

    const int incx_src   =src.grid_.inc(0);
    const int incy_src   =src.grid_.inc(1);
    const short nghosts_src=src.ghost_pt();

    int istart=0;
    int jstart=0;
    int kstart=0;
    if(dis=='g'){ // src is "global"
        istart=mype_env().my_mpi(0)*dim(0);
        jstart=mype_env().my_mpi(1)*dim(1);
        kstart=mype_env().my_mpi(2)*dim(2);
    }
 
    memset(uu_,0,grid_.sizeg()*sizeof(T));
    
    const T* const u_src=src.uu_
                             +kstart
                             +jstart*incy_src
                             +istart*incx_src;

    const size_t sdim2=dim(2)*sizeof(T);

    //if( dis=='d' )cout<<"expensive assign... "<<endl;
    for(int ix = 0;ix < dim_[0];ix++)
    {
        const int ix_dst = (ix+nghosts_dst)*incx_dst+nghosts_dst*(1+incy_dst);
        const int ix_src = (ix+nghosts_src)*incx_src+nghosts_src*(1+incy_src);

        for(int iy = 0;iy < dim_[1];iy++)
        {
            const int iy_dst = ix_dst+iy*incy_dst;
            const int iy_src = ix_src+iy*incy_src;
 
            memcpy(uu_+iy_dst, u_src+iy_src, sdim2);
        }
    }
    
    }
      
    updated_boundaries_=false;
}
template class GridFunc<double>;
template class GridFunc<float>;
} // namespace pb
