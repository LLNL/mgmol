// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SpreadPenaltyVolume.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"

using namespace std;

template <class T>
void SpreadPenaltyVolume<T>::addResidual(T& phi, T& res)
{
    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    assert(spreadf_ != nullptr);
    assert(spread_target_ >= 0.);
    assert(dampingFactor_ > 0.);

    const double total_volume_target = volume_target_ * phi.numst();

    // compute data to be used to evaluate spreads and centers
    spreadf_->computePositionMatrix(phi);

    // now compute spreads**2 and centers
    vector<Vector3D> centers;
    spreadf_->computeCenters(centers);

    // get spreads of overlapping orbitals
    vector<float> spread2;
    spreadf_->computeSpreads2(spread2);

    vector<float> local_spread2;
    spreadf_->computeLocalSpreads2(local_spread2);

    float sum_sigma6 = 0.;
    float sum_sigma3 = 0.;
    for (vector<float>::const_iterator it = local_spread2.begin();
         it != local_spread2.end(); ++it)
    {
        const float sigma3 = (*it) * sqrt(*it);

        sum_sigma3 += sigma3;
        sum_sigma6 += sigma3 * sigma3;
    }
    float sigmas[2] = { sum_sigma3, sum_sigma6 };
    mmpi.allreduce(&sigmas[0], 2, MPI_SUM);
    sum_sigma3 = sigmas[0];
    sum_sigma6 = sigmas[1];

    // if( onpe0 && ct.verbose>1 )cout<<"Average Spread
    // Volume="<<setprecision(4)<<sum_sigma3/phi.numst()<<endl;

    float coeff = 1. / (36. * sum_sigma6);
    coeff *= dampingFactor_;

    // add factor 3.*alpha*(sum_sigma3-total_volume_target) into coeff
    coeff *= 3. * alpha_ * (sum_sigma3 - total_volume_target);

    // get gids of overlapping orbitals
    const vector<int>& gids(spreadf_->getGids());

    // setup factors for corrections proportional to gradient of penalty spread
    // functional
    float max_spread2 = 0.;
    float min_spread2 = 1.e10;
    // float second_max_spread2=0.;
    // int gid_max_spread = -1;
    // int gid_min_spread = -1;
    // int   gid_second_max_spread=-1;

    vector<float> sigma;
    vector<int>::const_iterator gid_it = gids.begin();
    for (vector<float>::const_iterator it = spread2.begin();
         it != spread2.end(); ++it)
    {
        if (*it > max_spread2)
        {
            // second_max_spread2=max_spread2;
            max_spread2 = *it;
            // gid_second_max_spread=gid_max_spread;
            //       gid_max_spread = *gid_it;
        }
        if (*it < min_spread2)
        {
            // second_max_spread2=max_spread2;
            min_spread2 = *it;
            // gid_second_max_spread=gid_max_spread;
            // gid_min_spread = *gid_it;
        }

        float lm = sqrt(*it);
        sigma.push_back(lm);

        gid_it++;
    }

    // int rank;
    // mmpi.getRankMaxVal(max_spread2, &rank);
    // if(mype==rank)
    //{
    //    cout<<"PE "<<rank<<", Max. spread="<<sqrt(gmax_spread2)<<" for
    //    function "<<gid_max_spread<<endl;
    //}

    mmpi.allreduce(&max_spread2, 1, MPI_MAX);
    mmpi.allreduce(&min_spread2, 1, MPI_MIN);
    if (onpe0 && ct.verbose > 1)
        cout << "Max. spread = " << setprecision(4) << sqrt(max_spread2)
             << ", Min. spread = " << setprecision(4) << sqrt(min_spread2)
             << ", Average spread Volume = " << setprecision(4)
             << sum_sigma3 / phi.numst()
             << ", Target Volume = " << setprecision(4) << volume_target_
             << endl;

    // max_spread2=second_max_spread2;
    // gid_max_spread=gid_second_max_spread;

    // mmpi.getRankMaxVal(max_spread2, &rank);
    // mmpi.allreduce(&max_spread2,1,MPI_MAX);
    // if(mype==rank)
    //{
    //    cout<<"PE "<<rank<<", 2nd Max. spread="<<sqrt(max_spread2)<<" for
    //    function "<<gid_max_spread<<endl;
    //}

    computeAndAddResidualSpreadPenalty(sigma, coeff, centers, gids, phi, res);
}

// add to current residual "res" component proportional to gradient of penalty
// spread functional
template <class T>
void SpreadPenaltyVolume<T>::computeAndAddResidualSpreadPenalty(
    const vector<float>& sigma, const float eta,
    const vector<Vector3D>& centers, const vector<int>& gids, T& orbitals,
    T& res)
{
    assert(gids.size() == centers.size());

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    const int subdivx      = mymesh->subdivx();

    vector<float> lagrangemult;
    for (vector<float>::const_iterator it = sigma.begin(); it != sigma.end();
         ++it)
    {
        float lm = (*it) * (*it) * (*it);
        lagrangemult.push_back(lm);
    }

    // build maps: gid -> data
    map<int, Vector3D> gids2centers;
    short i = 0;
    for (vector<int>::const_iterator it = gids.begin(); it != gids.end(); ++it)
    {
        gids2centers.insert(pair<int, Vector3D>(*it, centers[i]));
        i++;
    }

    // build maps: gid -> data
    map<int, float> gids2sigma;
    i = 0;
    for (vector<int>::const_iterator it = gids.begin(); it != gids.end(); ++it)
    {
        gids2sigma.insert(pair<int, float>(*it, sigma[i]));
        i++;
    }

    map<int, float> gids2lagrangemult;
    i = 0;
    for (vector<int>::const_iterator it = gids.begin(); it != gids.end(); ++it)
    {
        gids2lagrangemult.insert(pair<int, float>(*it, lagrangemult[i]));
        i++;
    }

    const unsigned dim[3] = { mygrid.dim(0), mygrid.dim(1), mygrid.dim(2) };
    const int loc_length  = dim[0] / subdivx;
    const int incx        = dim[1] * dim[2];
    const int incy        = dim[2];

    const float inv_2pi = 0.5 * M_1_PI;
    float alphax        = mygrid.ll(0) * inv_2pi;
    float alphay        = mygrid.ll(1) * inv_2pi;
    float alphaz        = mygrid.ll(2) * inv_2pi;

    vector<float> sinx;
    vector<float> siny;
    vector<float> sinz;
    vector<float> cosx;
    vector<float> cosy;
    vector<float> cosz;
    mygrid.getSinCosFunctions(sinx, siny, sinz, cosx, cosy, cosz);

    const float coeffx = 2. * M_PI / mygrid.ll(0);
    const float coeffy = 2. * M_PI / mygrid.ll(1);
    const float coeffz = 2. * M_PI / mygrid.ll(2);

    const float pbound = 0.5;
    const float mbound = -1. * pbound;

    const vector<vector<int>>& global_indexes(orbitals.getOverlappingGids());

    for (short icolor = 0; icolor < orbitals.chromatic_number(); icolor++)
    {
        const ORBDTYPE* const iphi = orbitals.getPsi(icolor);
        ORBDTYPE* ires             = res.getPsi(icolor);

        for (short iloc = 0; iloc < subdivx; iloc++)
        {
            const int gid = global_indexes[iloc][icolor];

            if (gid > -1)
            {
                const float lambda = gids2lagrangemult[gid];
                const float sigma  = gids2sigma[gid];
                const Vector3D& center_gid(gids2centers[gid]);

                if (eta > 0.)
                    for (int ix = loc_length * iloc;
                         ix < loc_length * (iloc + 1); ix++)
                    {
                        const float sx0 = sin(coeffx * center_gid[0]) * alphax;
                        const float cx0 = cos(coeffx * center_gid[0]) * alphax;

                        const float x2 = (cosx[ix] - cx0) * (cosx[ix] - cx0)
                                         + (sinx[ix] - sx0) * (sinx[ix] - sx0);

                        for (unsigned int iy = 0; iy < dim[1]; iy++)
                        {
                            const float sy0
                                = sin(coeffy * center_gid[1]) * alphay;
                            const float cy0
                                = cos(coeffy * center_gid[1]) * alphay;

                            const float y2
                                = (cosy[iy] - cy0) * (cosy[iy] - cy0)
                                  + (siny[iy] - sy0) * (siny[iy] - sy0);

                            for (unsigned int iz = 0; iz < dim[2]; iz++)
                            {
                                const int index = ix * incx + iy * incy + iz;
                                const float cz0
                                    = cos(coeffz * center_gid[2]) * alphaz;
                                const float sz0
                                    = sin(coeffz * center_gid[2]) * alphaz;

                                const float z2
                                    = (cosz[iz] - cz0) * (cosz[iz] - cz0)
                                      + (sinz[iz] - sz0) * (sinz[iz] - sz0);

                                float val
                                    = eta
                                      * (2. * sigma * (x2 + y2 + z2) - lambda);
                                // if( fabs(val)>1. )cout<<"val="<<val<<endl;
                                val = val > pbound ? pbound : val;
                                val = val < mbound ? mbound : val;

                                ires[index] -= val * iphi[index];
                            }
                        }
                    }
            }
        }
    }
}

template <class T>
double SpreadPenaltyVolume<T>::evaluateEnergy(const T& phi)
{
    assert(spreadf_ != nullptr);
    assert(spread_target_ >= 0.);
    assert(dampingFactor_ > 0.);

    // compute data to be used to evaluate spreads and centers
    spreadf_->computePositionMatrix(phi);

    // now compute centers
    vector<Vector3D> centers;
    spreadf_->computeCenters(centers);

    // get spreads of functions centered in subdomain
    vector<float> spread2;
    spreadf_->computeLocalSpreads2(spread2);

    double sum_sigma3 = 0.;
    for (vector<float>::const_iterator it = spread2.begin();
         it != spread2.end(); ++it)
    {
        sum_sigma3 += ((*it) * sqrt(*it));
    }
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&sum_sigma3, 1, MPI_SUM);

    double delta = sum_sigma3 - volume_target_;
    if (delta > 0)
        return alpha_ * delta * delta;
    else
        return 0.;
}

template class SpreadPenaltyVolume<LocGridOrbitals<ORBDTYPE>>;
template class SpreadPenaltyVolume<ExtendedGridOrbitals<ORBDTYPE>>;
