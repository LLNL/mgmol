// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Grid.h,v 1.11 2010/01/28 22:56:47 jeanluc Exp $
#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>

#include "PEenv.h"
#include "Vector3D.h"
#include "../Control.h"

namespace pb{

class Grid
{
private:
    const PEenv&  mype_env_;

    int     dim_[3]; // local dim (on PE)
    int     gdim_[3];// global dim
    int     size_;
    int     sizeg_;
    long    gsize_;
    short   ghost_pt_;
    double  hgrid_[3];
    double  vel_;       // volume element
    int     inc_[3];
    double  ll_[3];
    double  start_[3];
    int     istart_[3];
    double  origin_[3];
    short   level_;
    
    bool    active_;
    
public:
    
    Grid(const double origin[3], 
         const double lattice[3], 
         const int ngpts[3], 
	 const PEenv& mype_env,
	 const short nghosts=1, const short level=0);
         
    // copy constructor
    Grid(const Grid&, const short nghosts=-1);

    Grid& operator=(const Grid&);

    bool active() const{ return active_;}
    const PEenv& mype_env() const{ return mype_env_; }
    int dim(const short i) const
    { 
        assert( dim_[i]>0 );
        assert( dim_[i]<10000 );
        return dim_[i]; 
    }
    int gdim(const short i) const{ return gdim_[i]; }
    int inc(const short i) const{ return inc_[i]; }
    double ll(const short i) const{ return ll_[i]; }
    double origin(const short i) const{ return origin_[i]; }
    double start(const short i) const{ return start_[i]; }
    int istart(const short i) const{ return istart_[i]; }

    int size() const{ return size_; }
    int sizeg() const{ return sizeg_; }
    long gsize() const{ return gsize_; }
    short ghost_pt() const
    { 
        assert( ghost_pt_<10 );
        return ghost_pt_;
    }
    
    double hgrid(short i) const{ return hgrid_[i]; }

    double vel() const{ return vel_; }


    short level() const{ return level_; }
    
    double hmax()const{
    	double hmax = hgrid_[0];
    	if(hgrid_[1] > hmax) hmax = hgrid_[1];
    	if(hgrid_[2] > hmax) hmax = hgrid_[2];
	return hmax;
    }
    
    double hmin()const{
    	double hmin = hgrid_[0];
    	if(hgrid_[1] < hmin) hmin = hgrid_[1];
    	if(hgrid_[2] < hmin) hmin = hgrid_[2];
	return hmin;
    }
    
    double anisotropy()const{
        return hmax()/hmin();
    }
    
    const Grid coarse_grid()const;
    const Grid replicated_grid(const PEenv&)const;
    
    Vector3D closestGridPt(Vector3D coords)const;

    ~Grid(){ 
            //cout<<"destructor for grid of dim"<<dim(0)<<"\n";
            }
    
    
    template <typename T>
    void getSinCosFunctions(
        std::vector<T>& sinx,
        std::vector<T>& siny,
        std::vector<T>& sinz,
        std::vector<T>& cosx,
        std::vector<T>& cosy,
        std::vector<T>& cosz)const;
};

} // namespace pb

#endif