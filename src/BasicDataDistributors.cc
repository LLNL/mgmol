// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/** Basic data distribution object calls */
#include "global.h"

#include "BasicDataDistributors.h"

#include "Control.h"
#include "LocalizationRegions.h"

DataDistribution* BasicDataDistributors::orbitalsProdWithHDistributor_
    = nullptr;
DataDistribution* BasicDataDistributors::gramMatrixDistributor_ = nullptr;
DataDistribution* BasicDataDistributors::centeredOrbitalsOverlapDistributor_
    = nullptr;

/* constructor */
BasicDataDistributors::BasicDataDistributors()
{
    orbitalsProdWithHDistributor_       = nullptr;
    gramMatrixDistributor_              = nullptr;
    centeredOrbitalsOverlapDistributor_ = nullptr;
}
/* destructor */
BasicDataDistributors::~BasicDataDistributors()
{
    delete orbitalsProdWithHDistributor_;
    orbitalsProdWithHDistributor_ = nullptr;
    delete gramMatrixDistributor_;
    gramMatrixDistributor_ = nullptr;
    delete centeredOrbitalsOverlapDistributor_;
    centeredOrbitalsOverlapDistributor_ = nullptr;
}
/* initialize */
void BasicDataDistributors::initialize(
    LocalizationRegions* lrs, const pb::PEenv& myPEenv, const double domain[])
{
    Control& ct(*(Control::instance()));
    // define extents for data distribution objects
    double spread_sH             = 3. * (*lrs).max_radii();
    double inverse_spread_radius = (*lrs).max_radii();

    assert(gramMatrixDistributor_ == nullptr);
    assert(centeredOrbitalsOverlapDistributor_ == nullptr);
    assert(orbitalsProdWithHDistributor_ == nullptr);

    gramMatrixDistributor_
        = new DataDistribution("gram", ct.spread_radius, myPEenv, domain);
    centeredOrbitalsOverlapDistributor_
        = new DataDistribution("invS", inverse_spread_radius, myPEenv, domain);
    orbitalsProdWithHDistributor_
        = new DataDistribution("prodH", spread_sH, myPEenv, domain);
}
