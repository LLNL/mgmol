// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef SUBCELL_H
#define SUBCELL_H

#include "Grid.h"
#include "Vector3D.h"

class SubCell
{
    double   inner_radius_;
    double   outer_radius_;
    std::vector<Vector3D> centers_;
    const Vector3D   cell_dimensions_;
    Vector3D   subcell_dimensions_;
    short subdivx_;
    std::vector<std::vector<double> > ll_;
    std::vector<std::vector<double> > ur_;

public:
    SubCell(const pb::Grid& fine_grid, const short subdivx, const short level):
        cell_dimensions_(fine_grid.ll(0),fine_grid.ll(1),fine_grid.ll(2)),
        subdivx_(subdivx)
    {
        assert( subdivx>0 );
        assert( level>=0 );
        
        double hgrid[3]={fine_grid.hgrid(0),fine_grid.hgrid(1),fine_grid.hgrid(2)};
        short l=0;
        while(l<level){
            for(short i=0;i<3;i++)
                hgrid[i]*=2.;
            l++;
        }

        const double start[3]={fine_grid.start(0),
                               fine_grid.start(1),
                               fine_grid.start(2)};
    
        const int dim[3]={(fine_grid.dim(0)>>level)/subdivx,
                          (fine_grid.dim(1)>>level),
                          (fine_grid.dim(2)>>level)};
        
        assert( dim[0]<=fine_grid.dim(0) );
        assert( dim[1]<=fine_grid.dim(1) );
        assert( dim[2]<=fine_grid.dim(2) );
        
        subcell_dimensions_[0]=(dim[0]-1)*hgrid[0];
        subcell_dimensions_[1]=(dim[1]-1)*hgrid[1],
        subcell_dimensions_[2]=(dim[2]-1)*hgrid[2];
        
        // subcell inner_radius
        // distance between subcell center and corners
        outer_radius_ = sqrt( 0.25*( subcell_dimensions_[0]*subcell_dimensions_[0]
                                    +subcell_dimensions_[1]*subcell_dimensions_[1] 
                                    +subcell_dimensions_[2]*subcell_dimensions_[2] ) );
        inner_radius_=std::min(subcell_dimensions_[0],subcell_dimensions_[1]);
        inner_radius_=std::min(inner_radius_,subcell_dimensions_[2]);
        inner_radius_*=0.5;
        
        centers_.resize(subdivx_);
        ll_.resize(subdivx_);
        ur_.resize(subdivx_);
        for(short iloc=0;iloc<subdivx_;iloc++){
            centers_[iloc][0] = start[0]+0.5*subcell_dimensions_[0]+iloc*dim[0]*hgrid[0];
            centers_[iloc][1] = start[1]+0.5*subcell_dimensions_[1];
            centers_[iloc][2] = start[2]+0.5*subcell_dimensions_[2];
            ll_[iloc].resize(3);
            ll_[iloc][0]=start[0]+iloc*dim[0]*hgrid[0];
            ll_[iloc][1]=start[1];
            ll_[iloc][2]=start[2];
            ur_[iloc].resize(3);
            for(short i=0;i<3;i++)
                ur_[iloc][i]=ll_[iloc][i]+subcell_dimensions_[i];
        }
    }
    
    double outerRadius()const
    {
        return outer_radius_;
    }
    
    // does a sphere of radius "radius" centered at "center" overlap with the subdomain?
    bool spherePossibleOverlap(const Vector3D& center, const double radius, const short iloc,
                               const short bcPoisson[3])const
    {
        // first find point in cell closest to center
        Vector3D vdcenter = centers_[iloc].vminimage(center, cell_dimensions_, bcPoisson);
        Vector3D vcorner;
        for(short i=0;i<3;i++)
        {
            vcorner[i] = vdcenter[i]<0. ? std::max(-0.5*subcell_dimensions_[i],vdcenter[i])
                                        : std::min( 0.5*subcell_dimensions_[i],vdcenter[i]);
        }

        // then check if this point is within radius of center
        for(short i=0;i<3;i++)
        {
            vdcenter[i]-=vcorner[i];
        }
        const double dcenter = length(vdcenter);

        return ( dcenter <= radius );
    }
    
    double distance(const Vector3D& center, const short iloc,
                    const short bcPoisson[3])const
    {
        return centers_[iloc].minimage(center, cell_dimensions_, bcPoisson);
    }
   
    bool includes(const Vector3D& point, const short iloc)const
    {
        const double tol=1.e-8;
        bool result=true;
        for(short i=0;i<3;i++)
        {
            result = result && (point[i]>=ll_[iloc][i]-tol);
            result = result && (point[i]<=ur_[iloc][i]+tol);
        }
        return result;
    }
 
    void printCenter(const short iloc, std::ostream& os)
    {
        os<<centers_[iloc];
        os<<std::endl;
    }
    
    void print(const short iloc, std::ostream& os)
    {
        os<<"Subdomain ("<<ll_[iloc][0]<<","<<ll_[iloc][1]<<","<<ll_[iloc][2]<<")-("
                         <<ur_[iloc][0]<<","<<ur_[iloc][1]<<","<<ur_[iloc][2]<<")"
                         <<std::endl;
    }
};

#endif
