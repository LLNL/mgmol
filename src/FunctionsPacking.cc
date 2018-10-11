// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "FunctionsPacking.h"
#include "Index.h"
#include "MPIdata.h"
#include "LocalizationRegions.h"
#include "SymmetricMatrix.h"
#include "coloring.h"

#include <list>
#include <cassert>
#include <numeric>
#include <cstdio>
#include <cstring>
using namespace std;


FunctionsPacking::FunctionsPacking(const FunctionsPacking& p)
{
    centers_    =p.centers_;
    gid2color_      =p.gid2color_;
    global_size_       =p.global_size_;
    num_colors_ =p.num_colors_;    
    loc_regions_=p.loc_regions_;
}

// Setup based on
// Recursive Largest First (RLF) algorithm
void FunctionsPacking::setup()
{
    Control& ct = *(Control::instance());
    list< list<int> > colored_gids;

    std::vector<int> gids;
    if( ct.globalColoring() )
        lrs_->getGidsGlobal(gids);
    else
        gids=lrs_->getOverlapGids();

    //SymmetricMatrix<char> orbi_overlap(numst_,comm_);
    SymmetricMatrix<char>* orbi_overlap;
    if( ct.RLFColoring() )
        orbi_overlap=new SymmetricMatrix<char>(gids,comm_);
    else
        orbi_overlap=new SymmetricMatrix<char>(gids.size(),comm_);
    if( ct.globalColoring() ){
        if( onpe0 && ct.verbose>1 )
            (*MPIdata::sout)<<" PACK STATES: Global coloring..."<<endl;
        initOrbiOverlapGlobal(0, *orbi_overlap);
    }else{
        if( onpe0 && ct.verbose>1 )
            (*MPIdata::sout)<<" PACK STATES: Local coloring..."<<endl;
        initOrbiOverlapLocal(0, *orbi_overlap);
    }

    global_size_=orbi_overlap->dimension();

    getColors(*orbi_overlap, colored_gids);

    setLocRegions(*lrs_,colored_gids);

    delete orbi_overlap;
}

//compute map "gid2color_" that maps gids to slots
void FunctionsPacking::getColors(const SymmetricMatrix<char>& overlaps, 
                                 list< list<int> >& colored_gids)
{
    assert( overlaps.dimension()<100000 );
    Control& ct = *(Control::instance());
    
    if( onpe0 && ct.verbose>0 )
    {
        (*MPIdata::sout)<<"setup FunctionsPacking for size="
                        <<global_size_<<endl;
    }

    if( !ct.RLFColoring() )
    {
        greedyColor(overlaps, colored_gids,
                    (ct.verbose>0 && onpe0), (*MPIdata::sout) );
    }else{
        colorRLF(overlaps, colored_gids,
                 (ct.verbose>0 && onpe0), (*MPIdata::sout) );
    }

    num_colors_=(short)colored_gids.size();
    if( onpe0 && ct.verbose>0 )
        (*MPIdata::sout)<<"FunctionsPacking::num_colors_= "
                        <<num_colors_<<endl; 

    assert( num_colors_<10000 );
   
    vector<short> nb_orb(num_colors_);

    // this is where we actually attribute gids to colors to slots
    short color=0;

    // loop over colored_gids
    list<list<int> >::const_iterator ic=colored_gids.begin();
    while( ic!=colored_gids.end() ){
        nb_orb[color]=(short)ic->size();
        //(*MPIdata::sout)<<"nb_orb["<<color<<"]="<<nb_orb[color]<<endl;

        // loop over functions of same color
        list<int>::const_iterator ii=ic->begin();
        while( ii!=ic->end() ){
            int gid=*ii; 
            gid2color_.insert( pair<int,short>(gid,color) );
            //if( onpe0 )
            //    (*MPIdata::sout)<<"Orbital "<<(*ii)<<" of color "<<gid2color_[(*ii)]<<endl;

            ii++;
        }
        color++;
        ic++;
    }
    
    // final check
    ic=colored_gids.begin();
    while( ic!=colored_gids.end() ){
        list<int>::const_iterator ii=ic->begin();
        while( ii!=ic->end() ){
            assert( gid2color_.find(*ii)->second<num_colors_ );
            ii++;
        }
        ic++;
    }

    int ntot_orb=accumulate(nb_orb.begin(),nb_orb.end(),0);
    
    if( ntot_orb!=global_size_ ){
        if( onpe0 )(*MPIdata::sout)<<ntot_orb<<" orbitals instead of "
                                   <<global_size_<<endl;
        exit(0);
    }

}

int FunctionsPacking::getAllCentersAndRadii4color(const short color, vector<double>& centers_and_radii)const
{    
    vector<double> local_centers_and_radii;
    getLocCentersAndRadii4color(color,local_centers_and_radii);
    
    centers_and_radii.clear();
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allGatherV(local_centers_and_radii,centers_and_radii);
    
    return centers_and_radii.size()/4;
}

void FunctionsPacking::setLocRegions(const LocalizationRegions& lrs,
                                     const list< list<int> >& colored_gids )
{
    loc_regions_.clear();
    Control& ct = *(Control::instance());
    
    if( ct.verbose>0 )
        printWithTimeStamp(" FunctionsPacking::setLocRegions",cout);
    
    list<list<int> >::const_iterator ic=colored_gids.begin();
    short color=0;
    if( lrs.globalNumLRs()>1 ){
        // loop over colored_gids
        while( ic!=colored_gids.end() ){

            // loop over functions of same color
            list<int>::const_iterator ii=ic->begin();
            while( ii!=ic->end() ){
                const int gid=*ii;
            
                assert( gid>=0 );
                assert( gid<(int)lrs.globalNumLRs() );
                if( lrs.radius(gid)>0.1 ) // radius=-1 if gid not know locally
                {
                    Sphere lr;
                    lr.center=lrs.getCenter(gid);
                    lr.radii =lrs.radius(gid);
                    assert( lr.radii>0. );
                    lr.gid   =gid;
                    loc_regions_.insert(pair<short,Sphere>(color,lr));
                }
                ii++;
            }
            color++;
            ic++;
        }
    }else{
        assert( lrs.globalNumLRs()==0 || lrs.globalNumLRs()==1 );
    
        // loop over colored_gids
        while( ic!=colored_gids.end() ){

            Sphere lr;
            lr.center=lrs.getCenter(0);
            lr.radii =lrs.radius(0);
            lr.gid   =*(ic->begin());
            loc_regions_.insert(pair<short,Sphere>(color,lr));
            
            color++;
            ic++;
        }
    }
}

// get localization centers and radii stored in function color
int FunctionsPacking::getLocCentersAndRadii4color(const short color,
                                                  vector<double>& data)const
{
    data.clear();
    multimap<short,Sphere>::const_iterator p  =loc_regions_.find(color);
    if( p==loc_regions_.end() )return 0;
    
    while( p!=loc_regions_.upper_bound(color) )
    {
        double arr[]={p->second.center[0],
                      p->second.center[1],
                      p->second.center[2],
                      p->second.radii};
        //data.push_back( p->second.center[0] );
        //data.push_back( p->second.center[1] );
        //data.push_back( p->second.center[2] );
        //data.push_back( p->second.radii );
        data.insert(data.end(),arr,arr+4);
        p++;
    }
    return data.size()/4;
}

// get localization centers and radii stored in function color
int FunctionsPacking::getLocGids4color(const short color, vector<int>& data)const
{
    data.clear();
    multimap<short,Sphere>::const_iterator p  =loc_regions_.find(color);
    if( p==loc_regions_.end() )return 0;

     while( p!=loc_regions_.upper_bound(color) )
     {
        data.push_back( p->second.gid );
        p++;
    }
    return data.size();
}

// get localization centers and radii stored in function color
int FunctionsPacking::getAllGids4color(const short color, vector<int>& gids)const
{
    vector<int> lgids;
    getLocGids4color(color,lgids);
    
    gids.clear();
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allGatherV(lgids,gids);
    
    return gids.size();
}

short FunctionsPacking::checkOverlapLRs(const int gid1, const int gid2)const
{
    return lrs_->overlap(gid1,gid2) ? 1 : 0;
}

//  Initialize the array orbi_overlap_[j][i] telling if 
//  the orbitals i and j overlap somewhere (on any PE or region)
// at level "level"
void FunctionsPacking::initOrbiOverlapLocal(const short level, SymmetricMatrix<char>& orbi_overlap)
{
    Control& ct = *(Control::instance());
    if( onpe0 && ct.verbose>2 )
        (*MPIdata::sout)<<"LocGridOrbitals::initOrbiOverlapLocal() for level "<<level<<endl;

    const vector<int>& local_gids(lrs_->getOverlapGids());
    
    assert( orbi_overlap.dimension()==local_gids.size() );
    
    // Initialization loop for orbi_overlap over all the states 
    for(vector<int>::const_iterator it1 =local_gids.begin();
                                    it1!=local_gids.end();
                                    it1++)
    {
        const int gid1=*it1;
        orbi_overlap.setVal(gid1,gid1,1);
        for(vector<int>::const_iterator it2 =local_gids.begin();
                                        it2!=it1;
                                        it2++)
        {
            const int gid2=*it2;
            short oo1=checkOverlapLRs(gid1, gid2);
            orbi_overlap.setVal(gid1,gid2,oo1);
        }
        
    }
}

void FunctionsPacking::initOrbiOverlapGlobal(const short level, SymmetricMatrix<char>& orbi_overlap)
{
    Control& ct = *(Control::instance());
    const int dim=orbi_overlap.dimension();
    if( onpe0 && ct.verbose>2 ){
        (*MPIdata::sout)<<"LocGridOrbitals::initOrbiOverlapGlobal() for level "<<level<<endl;
        (*MPIdata::sout)<<"LocGridOrbitals::initOrbiOverlapGlobal() for size "<<dim<<endl;
    }
    
    // Initialization loop for orbi_overlap over all the states 
    for(int gid1=0;gid1<dim;gid1++)
    {
        orbi_overlap.setVal(gid1,gid1,1);
        for(int gid2=0;gid2<gid1;gid2++)
        {
            short oo1=checkOverlapLRs(gid1, gid2);
            //short oo1=checkOverlap(gid1, gid2,level);
            //if( oo1!=oo2 )
            //{
            //    cout<<"mype="<<mype<<", gid1="<<gid1<<", gid2="<<gid2<<", oo1="<<oo1<<", oo2="<<oo2<<endl;
            //}
            //orbi_overlap_->setVal(gid1,gid2,checkOverlap(gid1, gid2,level));
            orbi_overlap.setVal(gid1,gid2,oo1);
        }
        
    }
#ifdef USE_MPI
    orbi_overlap.mpiAllOr();

#endif
}
