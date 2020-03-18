// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
 * Class to hold basic data distribution objects that are
 * likely to be most commonly used.
 */
#ifndef _BASICDATADISTRIBUTORS_H_
#define _BASICDATADISTRIBUTORS_H_

#include "DataDistribution.h"
class LocalizationRegions;

class BasicDataDistributors
{

    /* Data distribution objects */
    static DataDistribution* orbitalsProdWithHDistributor_;
    static DataDistribution* gramMatrixDistributor_;
    static DataDistribution* centeredOrbitalsOverlapDistributor_;

public:
    BasicDataDistributors();
    ~BasicDataDistributors();
    void initialize(LocalizationRegions* lrs, const pb::PEenv& myPEenv,
        const double domain[]);
    static DataDistribution* centeredOrbitalsOverlapDistributor()
    {
        return centeredOrbitalsOverlapDistributor_;
    }
    static DataDistribution* gramMatrixDistributor()
    {
        return gramMatrixDistributor_;
    }
    static DataDistribution* orbitalsProdWithHDistributor()
    {
        return orbitalsProdWithHDistributor_;
    }
};
#endif
