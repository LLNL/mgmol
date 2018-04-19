// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
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

void colorRLF(const SymmetricMatrix<char>& overlaps, 
           list< list<int> >& colored_gids)
{
    Control& ct = *(Control::instance());
    const int dim=overlaps.dimension();
    if( onpe0 && ct.verbose>0 )
    {
        (*MPIdata::sout)<<"Uses Recursive Largest First (RLF) algorithm"
                        <<" for problem of size "<<dim<<endl;
    }

    colored_gids.clear();
 
    // initialize list of "uncolored" gids
    list<int> uncolored;
    const vector<int>& gids(overlaps.gids());
    for(int i=0;i<dim;i++)uncolored.push_back(gids[i]);
    
    while( uncolored.size()>0 ){
    
        // compute degrees of vertices
        multimap<int,int> degrees;
        list<int>::const_iterator pj=uncolored.begin();
        while( pj!=uncolored.end() ){
            const int j0=(*pj);
            int degree=0;
            list<int>::const_iterator pi=uncolored.begin();
            while( pi!=uncolored.end() ){
                degree+=overlaps.val(*pi,j0);
                pi++;
            }
            if( degree )degrees.insert(pair<int,int>(degree,j0));
            pj++;
        }
            
        // find the 1st color: largest degree
        multimap<int,int>::const_iterator p0=degrees.end();
        p0--;        
        const int v0=p0->second;
        //if( onpe0 )
        //    (*MPIdata::sout)<<"Color "<<v0<<" of degree "<<p0->first<<endl;
        uncolored.remove(v0);
        
        // start list of functions of same color
        list<int> vi;
        vi.push_back(v0);
        
        // compute set u of vertices adjacents to v0
        // compute set v of vertices non-adjacents to v0
        list<int> u;
        list<int> v;
        pj=uncolored.begin();
        while( pj!=uncolored.end() ){
            int i0=(*pj);
            if( overlaps.val(i0,v0) ){
                u.push_back(i0);
            }else{
                v.push_back(i0);
            }
            pj++;
        }
        
        // try to color another vertex with color v0
        while( v.size()>0 ){
            
            // compute degrees in u for functions in v
            multimap<int,int> udegrees;
            list<int>::const_iterator pv=v.begin();
            while( pv!=v.end() ){
                const int j0=(*pv);
                int degree=0;
                list<int>::const_iterator pu=u.begin();
                while( pu!=u.end() ){
                    int i0=(*pu);
                    degree+=overlaps.val(i0,j0);
                    pu++;
                }
                udegrees.insert(pair<int,int>(degree,j0));
                pv++;
            }
                
            // find the largest degree in udegrees
            assert( udegrees.size()>0 );
            multimap<int,int>::const_iterator pu0=udegrees.end();
            pu0--;
            const int u0=pu0->second;
        
            // use same color for v0 and u0
            vi.push_back(u0);
            uncolored.remove(u0);
        
#ifndef NDEBUG
            if( onpe0 && ct.verbose>0 ){
                (*MPIdata::sout)<<"PE 0: Same color for orbitals "
                                <<*vi.begin()<<" and "<<u0<<endl;
                assert( !overlaps.val(*vi.begin(),u0) );
            }
#endif        
        
            // add to list u adjacents to u0
            list<int> vremove;
            pj=v.begin();
            while( pj!=v.end() ){
                const int i0=(*pj);
                if( overlaps.val(i0,u0) ){
                    u.push_back(i0);
                    vremove.push_back(i0);
                }
                pj++;
            }
            
            // remove from list v adjacents to u0
            pj=vremove.begin();
            while( pj!=vremove.end() ){
                v.remove(*pj);
                pj++;
            }

        }
        
        colored_gids.push_back(vi);
    }     
}

/* 
   quicksort of the elements in a. elements in b are permuted 
   according to the resulting sorted a.
*/
void quicksortI_h2l (int *a, int *b, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
    int i=lo, j=hi;
    int v;
    int x=ceil(a[(lo+hi)/2]);
    int q;

    //  partition
    do
    {    
        while (a[i]>x) i++; 
        while (a[j]<x) j--;
        if (i<=j)
        {
            v=a[i]; a[i]=a[j]; a[j]=v;
            q=b[i]; b[i]=b[j]; b[j]=q;
            i++; j--;
        }
    } while (i<=j);

    //  recursion
    if (lo<j) quicksortI_h2l(a, b, lo, j);
    if (i<hi) quicksortI_h2l(a, b, i, hi);
}   

void greedyColor(const SymmetricMatrix<char>& overlaps, 
           list< list<int> >& colored_gids)
{
    if( onpe0 )
    {
        (*MPIdata::sout)<<"Uses greedy algorithm"<<endl;
    }

   const int dim = overlaps.dimension();
   vector<int>vecia;
   vector<int>vecja;
   vecia.reserve(dim+1);
   vecja.reserve(dim);

   int *nnzrow = new int[dim];
   int *iord = new int[dim];
   
   /* nonzero pattern of overlap matrix */
   overlaps.getNNZPattern(vecia, vecja, nnzrow);
   
   for(int i=0; i<dim; i++)
      iord[i] = i;
   
   /* sort degree of nodes from hi to lo */
   quicksortI_h2l(nnzrow, iord, 0, dim-1);
   /* perform greedy multicoloring */
   int num_colors=0;
   vector<int>::iterator ia = vecia.begin();
   vector<int>::iterator ja = vecja.begin();   
   /* overwrite nnzrow array -- no longer needed */
   int *colors = nnzrow;
   greedyMC(dim, &(*ja), &(*ia), &num_colors, iord, colors);
   
   /* populate list of colored_grids */
   for(int color=0; color<num_colors; color++)
   {
      list<int> llist;
      for(int i=0; i<dim; i++)
      {
         if(colors[i] == color)
            llist.push_back(i);
      }
      colored_gids.push_back(llist);
   }

   delete [] iord;
   delete [] colors;

   return;
}  
/*----------------------------------------------------------------*/

FunctionsPacking::FunctionsPacking(const FunctionsPacking& p)
{
    centers_    =p.centers_;
    where_      =p.where_;
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

void FunctionsPacking::getColors(const SymmetricMatrix<char>& overlaps, 
                                 list< list<int> >& colored_gids)
{
    assert( overlaps.dimension()<100000 );
    Control& ct = *(Control::instance());
    
    if( onpe0 && ct.verbose>0 )
    {
        (*MPIdata::sout)<<"setup FunctionsPacking for size="<<global_size_<<endl;
    }

    if( !ct.RLFColoring() )
    {
        greedyColor(overlaps, colored_gids);
    }else{
        colorRLF(overlaps, colored_gids);
    }

    num_colors_=(short)colored_gids.size();
    if( onpe0 && ct.verbose>0 )(*MPIdata::sout)<<" num_colors_= "<<num_colors_<<endl; 

    assert( num_colors_<10000 );
   
    vector<short> nb_orb(num_colors_);

    // attribute colors
    short color=0;

    // loop over colored_gids
    list<list<int> >::const_iterator ic=colored_gids.begin();
    while( ic!=colored_gids.end() ){
        nb_orb[color]=(short)ic->size();
        //(*MPIdata::sout)<<"nb_orb["<<color<<"]="<<nb_orb[color]<<endl;

        // loop over functions of same color
        list<int>::const_iterator ii=ic->begin();
        while( ii!=ic->end() ){
            
            where_.insert( pair<int,short>(*ii,color) );
            //if( onpe0 )
            //    (*MPIdata::sout)<<"Orbital "<<(*ii)<<" of color "<<where_[(*ii)]<<endl;

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
            assert( where_.find(*ii)->second<num_colors_ );
            ii++;
        }
        ic++;
    }

    int ntot_orb=accumulate(nb_orb.begin(),nb_orb.end(),0);
    
    if( ntot_orb!=global_size_ ){
        if( onpe0 )(*MPIdata::sout)<<ntot_orb<<" orbitals instead of "<<global_size_<<endl;
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
