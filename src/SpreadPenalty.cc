// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SpreadPenalty.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"
#include "memory_space.h"

#ifdef HAVE_MAGMA
using memory_space_type = MemorySpace::Device;
#else
using memory_space_type = MemorySpace::Host;
#endif

template <class T>
void SpreadPenalty<T>::addResidual(T& phi, T& res)
{
    Control& ct = *(Control::instance());

    assert(spreadf_ != nullptr);
    assert(spread2_target_ >= 0.);

    // compute data to be used to evaluate spreads and centers
    spreadf_->computePositionMatrix(phi);

    // compute centers (within computational domain)
    std::vector<Vector3D> centers;
    spreadf_->computeCenters(centers);

    // get spreads of overlapping orbitals
    std::vector<float> spread2;
    spreadf_->computeSpreads2(spread2);

    // get gids of overlapping orbitals
    const std::vector<int>& gids(spreadf_->getGids());

    // setup factors for corrections proportional to gradient
    // of penalty spread functional
    std::vector<float> factors;
    factors.reserve(spread2.size());
    float max_spread2 = 0.;
    // float second_max_spread2=0.;
    // int gid_max_spread = -1;
    // int   gid_second_max_spread=-1;

    for (auto it : spread2)
    {
        if (it > max_spread2)
        {
            max_spread2 = it;
        }

        // Compute appropriate spread penalty factor
        float coeff = computeSpreadPenaltyFactor(it);

        factors.push_back(coeff);
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    // int rank;
    // mmpi.getRankMaxVal(max_spread2, &rank);
    // if(mype==rank)
    //{
    //    cout<<"PE "<<rank<<", Max. spread="<<sqrt(gmax_spread2)<<" for
    //    function "<<gid_max_spread<<endl;
    //}

    mmpi.allreduce(&max_spread2, 1, MPI_MAX);
    if (onpe0 && ct.verbose > 1)
        std::cout << "Max. spread = " << std::setprecision(4)
                  << sqrt(max_spread2) << std::endl;

    // max_spread2=second_max_spread2;
    // gid_max_spread=gid_second_max_spread;

    // mmpi.getRankMaxVal(max_spread2, &rank);
    // mmpi.allreduce(&max_spread2,1,MPI_MAX);
    // if(mype==rank)
    //{
    //    cout<<"PE "<<rank<<", 2nd Max. spread="<<sqrt(max_spread2)<<" for
    //    function "<<gid_max_spread<<endl;
    //}

    computeAndAddResidualSpreadPenalty(
        spread2, factors, centers, gids, phi, res);
}

template <class T>
float SpreadPenalty<T>::computeSpreadPenaltyFactor(const float spread2)
{
    // compute 2*F[\phi]-sigma_0^2
    //(alpha ignored since it cancels out with the one in eta)
    float factor = (spread2 - spread2_target_);

    factor = factor > 0. ? (2. * factor) : 0.;

    factor /= (8. * spread2 * spread2);

    factor *= dampingFactor_;

    return factor;
}

template <class T>
float SpreadPenalty<T>::computeSpreadPenaltyFactorXLBOMD(const float spread2)
{
    float factor = (spread2 - spread2_target_);

    factor = factor > 0. ? (4. * alpha_ * factor) : 0.;

    // factor*=dampingFactor_;

    return factor;
}

// add to current residual "res" component proportional to gradient of penalty
// spread functional
template <class T>
void SpreadPenalty<T>::computeAndAddResidualSpreadPenalty(
    const std::vector<float>& lagrangemult, const std::vector<float>& factors,
    const std::vector<Vector3D>& centers, const std::vector<int>& gids,
    T& orbitals, T& res)
{
    assert(lagrangemult.size() == centers.size());
    assert(factors.size() == centers.size());
    assert(gids.size() == centers.size());

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    const int subdivx      = mymesh->subdivx();

    // build std::maps: gid -> data
    std::map<int, Vector3D> gids2centers;
    short i = 0;
    for (std::vector<int>::const_iterator it = gids.begin(); it != gids.end();
         ++it)
    {
        gids2centers.insert(std::pair<int, Vector3D>(*it, centers[i]));
        i++;
    }

    std::map<int, float> gids2lagrangemult;
    i = 0;
    for (std::vector<int>::const_iterator it = gids.begin(); it != gids.end();
         ++it)
    {
        gids2lagrangemult.insert(std::pair<int, float>(*it, lagrangemult[i]));
        i++;
    }

    std::map<int, float> gids2factors;
    i = 0;
    for (std::vector<int>::const_iterator it = gids.begin(); it != gids.end();
         ++it)
    {
        gids2factors.insert(std::pair<int, float>(*it, factors[i]));
        i++;
    }

    const unsigned dim[3] = { mygrid.dim(0), mygrid.dim(1), mygrid.dim(2) };
    const int loc_length  = dim[0] / subdivx;
    const int incx        = dim[1] * dim[2];
    const int incy        = dim[2];

    const float inv_2pi = 0.5 * M_1_PI;

    std::vector<float> sinx;
    std::vector<float> siny;
    std::vector<float> sinz;
    std::vector<float> cosx;
    std::vector<float> cosy;
    std::vector<float> cosz;
    mygrid.getSinCosFunctions(sinx, siny, sinz, cosx, cosy, cosz);

    const float l2rad[3] = { (float)(2. * M_PI / mygrid.ll(0)),
        (float)(2. * M_PI / mygrid.ll(1)), (float)(2. * M_PI / mygrid.ll(2)) };
    const float rad2l[3] = { (float)(mygrid.ll(0) * inv_2pi),
        (float)(mygrid.ll(1) * inv_2pi), (float)(mygrid.ll(2) * inv_2pi) };

    const float pbound = 0.5;
    const float mbound = -1. * pbound;

    const std::vector<std::vector<int>>& global_indexes(
        orbitals.getOverlappingGids());

    unsigned int const phi_size = mygrid.size();
    ORBDTYPE* phi_host_view
        = MemorySpace::Memory<ORBDTYPE, memory_space_type>::allocate_host_view(
            phi_size);
    ORBDTYPE* res_host_view
        = MemorySpace::Memory<ORBDTYPE, memory_space_type>::allocate_host_view(
            phi_size);

    for (short icolor = 0; icolor < orbitals.chromatic_number(); icolor++)
    {
        ORBDTYPE* iphi = orbitals.getPsi(icolor);
        ORBDTYPE* ires = res.getPsi(icolor);

        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            iphi, phi_size, phi_host_view);
        // initialize res_host_view with current residual as we are going
        // to add contributions to it
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            ires, phi_size, res_host_view);

        for (short iloc = 0; iloc < subdivx; iloc++)
        {
            const int gid = global_indexes[iloc][icolor];

            if (gid > -1)
            {
                const float lambda = gids2lagrangemult[gid];
                Vector3D center(gids2centers[gid]);
                const float factor = gids2factors[gid];

                // map center into [0,2pi]
                for (int i = 0; i < 3; i++)
                {
                    center[i] -= mygrid.origin(i);
                    center[i] *= l2rad[i];
                }
                // std::cout<<"center = "<<center<<std::endl;

                if (factor > 0.)
                    for (int ix = loc_length * iloc;
                         ix < loc_length * (iloc + 1); ix++)
                    {
                        const float sx0 = std::sin(center[0]) * rad2l[0];
                        const float cx0 = std::cos(center[0]) * rad2l[0];

                        const float x2 = (cosx[ix] - cx0) * (cosx[ix] - cx0)
                                         + (sinx[ix] - sx0) * (sinx[ix] - sx0);

                        for (unsigned int iy = 0; iy < dim[1]; iy++)
                        {
                            const float sy0 = std::sin(center[1]) * rad2l[1];
                            const float cy0 = std::cos(center[1]) * rad2l[1];

                            const float y2
                                = (cosy[iy] - cy0) * (cosy[iy] - cy0)
                                  + (siny[iy] - sy0) * (siny[iy] - sy0);

                            for (unsigned int iz = 0; iz < dim[2]; iz++)
                            {
                                const int index = ix * incx + iy * incy + iz;
                                const float cz0
                                    = std::cos(center[2]) * rad2l[2];
                                const float sz0
                                    = std::sin(center[2]) * rad2l[2];

                                const float z2
                                    = (cosz[iz] - cz0) * (cosz[iz] - cz0)
                                      + (sinz[iz] - sz0) * (sinz[iz] - sz0);

                                float val = factor * (x2 + y2 + z2 - lambda);
                                // if( fabs(val)>1. )cout<<"val="<<val<<endl;
                                val = val > pbound ? pbound : val;
                                val = val < mbound ? mbound : val;

                                res_host_view[index]
                                    -= val * phi_host_view[index];
                            }
                        }
                    }
            }
        }

        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_dev(
            res_host_view, phi_size, ires);
    }

    MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
        phi_host_view);
    MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
        res_host_view);
}

template <class T>
double SpreadPenalty<T>::evaluateEnergy(const T& phi)
{
    assert(spreadf_ != nullptr);
    assert(spread2_target_ >= 0.);
    assert(dampingFactor_ > 0.);

    // compute data to be used to evaluate spreads and centers
    spreadf_->computePositionMatrix(phi);

    // now compute centers
    std::vector<Vector3D> centers;
    spreadf_->computeCenters(centers);

    // get spreads of functions centered in subdomain
    std::vector<float> spread2;
    spreadf_->computeLocalSpreads2(spread2);

    double total_energy = 0.;

    for (std::vector<float>::const_iterator it = spread2.begin();
         it != spread2.end(); ++it)
    {
        double diff = (*it - spread2_target_);

        if (diff > 0.) total_energy += (diff * diff);
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&total_energy, 1, MPI_SUM);

    return alpha_ * total_energy;
}

template class SpreadPenalty<LocGridOrbitals<ORBDTYPE>>;
template class SpreadPenalty<ExtendedGridOrbitals<ORBDTYPE>>;
