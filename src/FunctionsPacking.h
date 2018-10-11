// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
// Class to allocate orbitals into "global" functions
#include "Vector3D.h"
#include "SymmetricMatrix.h"

#include<map>
#include<set>
#include<vector>
#include<list>

#ifndef FUNCTIONSPACKING_H
#define FUNCTIONSPACKING_H

class LocalizationRegions;

typedef struct Sphere
{
  Vector3D center;
  double   radii;
  int      gid;
} Sphere;


class FunctionsPacking
{
private:
    std::vector<Vector3D>              centers_;
    std::map<int,short>                gid2color_; // tells where one gid is allocated
    int                                global_size_;
    short                              num_colors_;
    
    LocalizationRegions*   lrs_;
    MPI_Comm    comm_;

    // orbitals centers (x,y,z) and radii for each global array
    // 1st index=color, 2nd index=center+radii+gid
    std::multimap<short,Sphere>  loc_regions_;
    
    void setLocRegions(const LocalizationRegions& lrs,
                       const std::list< std::list<int> >& colors );
    void getColors(const SymmetricMatrix<char>& overlaps, 
                   std::list< std::list<int> >& colors);
    short checkOverlapLRs(const int gid1, const int gid2)const;
    void initOrbiOverlapLocal(const short level, SymmetricMatrix<char>& orbi_overlap);
    void initOrbiOverlapGlobal(const short level, SymmetricMatrix<char>& orbi_overlap);

public:

    
    FunctionsPacking(LocalizationRegions* lrs, const MPI_Comm comm)
       :lrs_(lrs),
        comm_(comm)
    {
    };
    FunctionsPacking(const FunctionsPacking&);
    
    void setup();

    // return color of gid if exists locally, otherwise return -1
    short getColor(const int gid)const
    { 
        assert( gid>=0 );
        
        std::map<int,short>::const_iterator it = gid2color_.find(gid);
        if( it==gid2color_.end() )
           return -1;
        else
           return it->second;
    }
    void getPossibleColors(const Vector3D center, std::set<int>& colors)
    {        
        colors.clear();
        std::multimap<short,Sphere>::iterator p=loc_regions_.begin();
        while( p!=loc_regions_.end() )
        {
            if( p->second.center==center )colors.insert(p->first);
            p++;
        }
    }
    
    // get localization centers stored in function color
    int getLocCentersAndRadii4color(const short color,
        std::vector<double>& centers_and_radii)const;
    int getAllCentersAndRadii4color(const short color,
        std::vector<double>& centers_and_radii)const;

    int getLocGids4color(const short color, std::vector<int>& data)const;
    int getAllGids4color(const short color, std::vector<int>& data)const;
   

    int chromatic_number()const
    {
        return num_colors_;
    }

};

#endif
