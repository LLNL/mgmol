// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SpreadPenalty.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"

void SpreadPenalty::addResidual(LocGridOrbitals& phi,
                                LocGridOrbitals& res)
{
    addResidual(phi, res, false);
}

void SpreadPenalty::addResidual(LocGridOrbitals& phi,
                                LocGridOrbitals& res,
                                bool xlbomd)
{
    Control& ct = *(Control::instance());
    
    assert( spreadf_!=0 );
    assert( spread2_target_>=0. );
    
    
    // compute data to be used to evaluate spreads and centers
    spreadf_->computePositionMatrix(phi);
        
    // now compute spreads**2 and centers
    vector<Vector3D> centers;
    spreadf_->computeCenters(centers);
    
    // get spreads of overlapping orbitals
    vector<float> spread2;
    spreadf_->computeSpreads2(spread2);
    
    // get gids of overlapping orbitals
    const vector<int>& gids( spreadf_->getGids() );

    // setup factors for corrections proportional to gradient 
    // of penalty spread functional
    vector<float> factors;
    float max_spread2=0.;
    //float second_max_spread2=0.;
    int   gid_max_spread=-1;
    //int   gid_second_max_spread=-1;
    
    vector<int>::const_iterator gid_it=gids.begin();
    for(vector<float>::const_iterator it =spread2.begin();
                                      it!=spread2.end();
                                    ++it)
    {
        if( *it>max_spread2 )
        {
            //second_max_spread2=max_spread2;
            max_spread2=*it;
            //gid_second_max_spread=gid_max_spread;
            gid_max_spread=*gid_it;
        }
       
        // Compute appropriate spread penalty factor for either localized XLBOMD or penalized optimization.
        float coeff;
            coeff = computeSpreadPenaltyFactor(*it);
 
        factors.push_back( coeff );
        
        gid_it++;
    }
    
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        
    //int rank;
    //mmpi.getRankMaxVal(max_spread2, &rank);
    //if(mype==rank)
    //{
    //    cout<<"PE "<<rank<<", Max. spread="<<sqrt(gmax_spread2)<<" for function "<<gid_max_spread<<endl;
    //}
    
    mmpi.allreduce(&max_spread2,1,MPI_MAX);
    if( onpe0 && ct.verbose>1 )cout<<"Max. spread="<<setprecision(4)<<sqrt(max_spread2)<<endl;
    
    //max_spread2=second_max_spread2;
    //gid_max_spread=gid_second_max_spread;
    
    //mmpi.getRankMaxVal(max_spread2, &rank);    
    //mmpi.allreduce(&max_spread2,1,MPI_MAX);
    //if(mype==rank)
    //{
    //    cout<<"PE "<<rank<<", 2nd Max. spread="<<sqrt(max_spread2)<<" for function "<<gid_max_spread<<endl;
    //}

    computeAndAddResidualSpreadPenalty(spread2, factors, centers, gids, phi, res);
}

float SpreadPenalty::computeSpreadPenaltyFactor(const float spread2)
{
    //compute 2*F[\phi]-sigma_0^2
    //(alpha ignored since it cancels out with teh one in eta)
    float factor = (spread2 - spread2_target_);

    factor = factor>0. ? (2.*factor) : 0.;

    factor/=(8.*spread2*spread2);

    factor*=dampingFactor_;

    return factor;
}

float SpreadPenalty::computeSpreadPenaltyFactorXLBOMD(const float spread2)
{
    float factor = (spread2 - spread2_target_);

    factor = factor>0. ? (4.*alpha_*factor) : 0.;

    //factor*=dampingFactor_;

    return factor;
}

// add to current residual "res" component proportional to gradient of penalty spread functional
void SpreadPenalty::computeAndAddResidualSpreadPenalty(const vector<float>& lagrangemult,
                                                       const vector<float>& factors,
                                                       const vector<Vector3D>& centers,
                                                       const vector<int>& gids,
                                                       LocGridOrbitals& orbitals,
                                                       LocGridOrbitals& res)
{
    assert( lagrangemult.size()==centers.size() );
    assert( factors.size()==centers.size() );
    assert( gids.size()==centers.size() );

    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();  
    const int subdivx=mymesh->subdivx();
    
    // build maps: gid -> data
    map<int,Vector3D> gids2centers;
    short i=0;
    for(vector<int>::const_iterator it =gids.begin();
                                    it!=gids.end();
                                  ++it)
    {
        gids2centers.insert(pair<int,Vector3D>(*it,centers[i]));
        i++;
    }
    
    map<int,float> gids2lagrangemult;
    i=0;
    for(vector<int>::const_iterator it =gids.begin();
                                    it!=gids.end();
                                  ++it)
    {
        gids2lagrangemult.insert(pair<int,float>(*it,lagrangemult[i]));
        i++;
    }
    
    map<int,float> gids2factors;
    i=0;
    for(vector<int>::const_iterator it =gids.begin();
                                    it!=gids.end();
                                  ++it)
    {
        gids2factors.insert(pair<int,float>(*it,factors[i]));
        i++;
    }
    
    const int dim[3]={mygrid.dim(0),mygrid.dim(1),mygrid.dim(2)};
    const int loc_length = dim[0]/subdivx;
    const int incx=dim[1]*dim[2];
    const int incy=dim[2];

    const float inv_2pi=0.5*M_1_PI;
    float  alphax=mygrid.ll(0)*inv_2pi;
    float  alphay=mygrid.ll(1)*inv_2pi;
    float  alphaz=mygrid.ll(2)*inv_2pi;

    vector<float> sinx;
    vector<float> siny;
    vector<float> sinz;
    vector<float> cosx;
    vector<float> cosy;
    vector<float> cosz;
    mygrid.getSinCosFunctions(sinx,siny,sinz,cosx,cosy,cosz);
    
    const float coeffx=2.*M_PI/mygrid.ll(0);
    const float coeffy=2.*M_PI/mygrid.ll(1);
    const float coeffz=2.*M_PI/mygrid.ll(2);
    
    const float pbound=0.5;
    const float mbound=-1.*pbound;
    
    const vector<vector<int> >& global_indexes(orbitals.getGlobalIndexes());
    
    for(short icolor=0;icolor<orbitals.chromatic_number();icolor++)
    {
        const ORBDTYPE* const iphi=orbitals.getPhi(icolor);
        ORBDTYPE* ires=res.getPhi(icolor);

        for(short iloc=0;iloc<subdivx;iloc++)
        {
            const int gid=global_indexes[iloc][icolor];
            
            if(gid>-1)
            {
                const float lambda=gids2lagrangemult[gid];
                const Vector3D& center_gid( gids2centers[gid] );
                const float factor=gids2factors[gid];
                
                if( factor>0. )
                for(int ix = loc_length*iloc;ix < loc_length*(iloc+1);ix++)
                {
                    const float sx0=sin(coeffx*center_gid[0])*alphax;
                    const float cx0=cos(coeffx*center_gid[0])*alphax;
                    
                    const float x2=(cosx[ix]-cx0)*(cosx[ix]-cx0)+(sinx[ix]-sx0)*(sinx[ix]-sx0);
                    
                    for(int iy = 0;iy < dim[1];iy++)
                    {
                        const float sy0=sin(coeffy*center_gid[1])*alphay;
                        const float cy0=cos(coeffy*center_gid[1])*alphay;
                        
                        const float y2=(cosy[iy]-cy0)*(cosy[iy]-cy0)+(siny[iy]-sy0)*(siny[iy]-sy0);
                        
                        for(int iz = 0;iz < dim[2];iz++)
                        {
                            const int index=ix*incx+iy*incy+iz;
                            const float cz0=cos(coeffz*center_gid[2])*alphaz;
                            const float sz0=sin(coeffz*center_gid[2])*alphaz;
                            
                            const float z2=(cosz[iz]-cz0)*(cosz[iz]-cz0)+(sinz[iz]-sz0)*(sinz[iz]-sz0);
                            
                            float val=factor*(x2+y2+z2 -lambda );
                            //if( fabs(val)>1. )cout<<"val="<<val<<endl;
                            val = val > pbound ? pbound : val;
                            val = val < mbound ? mbound : val;
                            
                            ires[index] -= val*iphi[index];
                        
                        }
                    }
                }
            }
        }
    }

}

                                                   
double SpreadPenalty::evaluateEnergy(const LocGridOrbitals& phi)
{    
    assert( spreadf_!=0 );
    assert( spread2_target_>=0. );
    assert( dampingFactor_>0. );
    
    // compute data to be used to evaluate spreads and centers
    spreadf_->computePositionMatrix(phi);
        
    // now compute centers
    vector<Vector3D> centers;
    spreadf_->computeCenters(centers);
    
    // get spreads of functions centered in subdomain
    vector<float> spread2;
    spreadf_->computeLocalSpreads2(spread2);

    double total_energy=0.;
    
    for(vector<float>::const_iterator it =spread2.begin();
                                      it!=spread2.end();
                                    ++it)
    {
        double diff = (*it - spread2_target_);

        if( diff>0. )total_energy += (diff*diff);        
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&total_energy,1,MPI_SUM);
    
    return alpha_*total_energy;
}
