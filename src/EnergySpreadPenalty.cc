// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "EnergySpreadPenalty.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"

template <class T>
void EnergySpreadPenalty<T>::addResidual(T& phi, T& res)
{
    assert(spreadf_ != nullptr);
    assert(spread2_target_ >= 0.);

    // compute data to be used to evaluate spreads and centers
    spreadf_->computePositionMatrix(phi);

    // now compute spreads**2 and centers
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
    float max_spread2 = 0.;

    std::vector<int>::const_iterator gid_it = gids.begin();
    for (std::vector<float>::const_iterator it = spread2.begin();
         it != spread2.end(); ++it)
    {
        // compute 0.5*2*alpha*(F[\phi]-sigma_0^2)
        // factor 0.5 is to be consistent with residual of DFT
        // H*phi-...
        float coeff = (*it - spread2_target_);
        // if( coeff>0. )cout<<"spread = "<<sqrt(*it)<<endl;
        if (*it > max_spread2)
        {
            max_spread2 = *it;
        }
        coeff /= (8. * (*it) * (*it));

        coeff = coeff > 0. ? (coeff) : 0.;
        factors.push_back(coeff);

        gid_it++;
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    mmpi.allreduce(&max_spread2, 1, MPI_MAX);
    if (onpe0)
        std::cout << "Max. spread = " << std::setprecision(4)
                  << std::sqrt(max_spread2) << std::endl;

    computeAndAddResidualSpreadPenalty(
        spread2, factors, centers, gids, phi, res);
}

// add to current residual "res" component proportional to gradient of penalty
// spread functional factors equal 2*alpha(F[\phi]-sigma_0^2) (not included in
// lagrangemult) assumes orbitals are normalized
template <class T>
void EnergySpreadPenalty<T>::computeAndAddResidualSpreadPenalty(
    const std::vector<float>& spread2, const std::vector<float>& factors,
    const std::vector<Vector3D>& centers, const std::vector<int>& gids,
    T& orbitals, T& res)
{
    assert(spread2.size() == centers.size());
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

    std::map<int, float> gids2spread2;
    i = 0;
    for (std::vector<int>::const_iterator it = gids.begin(); it != gids.end();
         ++it)
    {
        gids2spread2.insert(std::pair<int, float>(*it, spread2[i]));
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
    float alphax        = mygrid.ll(0) * inv_2pi;
    float alphay        = mygrid.ll(1) * inv_2pi;
    float alphaz        = mygrid.ll(2) * inv_2pi;

    std::vector<float> sinx;
    std::vector<float> siny;
    std::vector<float> sinz;
    std::vector<float> cosx;
    std::vector<float> cosy;
    std::vector<float> cosz;
    mygrid.getSinCosFunctions(sinx, siny, sinz, cosx, cosy, cosz);

    const float coeffx = 2. * M_PI / mygrid.ll(0);
    const float coeffy = 2. * M_PI / mygrid.ll(1);
    const float coeffz = 2. * M_PI / mygrid.ll(2);

    const std::vector<std::vector<int>>& global_indexes(
        orbitals.getOverlappingGids());

    for (short icolor = 0; icolor < orbitals.chromatic_number(); icolor++)
    {
        const ORBDTYPE* const iphi = orbitals.getPsi(icolor);
        ORBDTYPE* ires             = res.getPsi(icolor);

        for (short iloc = 0; iloc < subdivx; iloc++)
        {
            const int gid = global_indexes[iloc][icolor];

            if (gid > -1)
            {
                const float fphi = gids2spread2[gid];
                const Vector3D& center_gid(gids2centers[gid]);
                const float factor = gids2factors[gid];

                if (factor > 0.)
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

                                float val = factor * (x2 + y2 + z2 - fphi);

                                ires[index] -= 2. * val * iphi[index];
                            }
                        }
                    }
            }
        }
    }
}

template <class T>
double EnergySpreadPenalty<T>::evaluateEnergy(const T& phi)
{
    assert(spreadf_ != nullptr);
    assert(spread2_target_ >= 0.);

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

template class EnergySpreadPenalty<LocGridOrbitals<ORBDTYPE>>;
template class EnergySpreadPenalty<ExtendedGridOrbitals<ORBDTYPE>>;
