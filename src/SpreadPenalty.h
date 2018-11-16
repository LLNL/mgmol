// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef SpreadPenalty_H
#define SpreadPenalty_H

class LocGridOrbitals;
#include "Control.h"
#include "SpreadPenaltyInterface.h"
#include "SpreadsAndCenters.h"

class SpreadPenalty : public SpreadPenaltyInterface
{
private:
    SpreadsAndCenters* spreadf_;

    double spread2_target_;
    double alpha_;
    double dampingFactor_;

    void computeAndAddResidualSpreadPenalty(const vector<float>& lagrangemult,
        const vector<float>& factors, const vector<Vector3D>& centers,
        const vector<int>& gids, LocGridOrbitals& orbitals,
        LocGridOrbitals& res);

    float computeSpreadPenaltyFactor(const float spread2);
    float computeSpreadPenaltyFactorXLBOMD(const float spread2);

public:
    SpreadPenalty(SpreadsAndCenters* spreadf, const double spread_target,
        const double alpha, const double dampingFactor)
        : spreadf_(spreadf),
          spread2_target_(spread_target * spread_target),
          alpha_(alpha),
          dampingFactor_(dampingFactor)
    {
        assert(spread2_target_ >= 0.);
        assert(alpha_ >= 0.);
    }

    // add penalty functional contribution to residual
    void addResidual(LocGridOrbitals& phi, LocGridOrbitals& res);

    void addResidual(LocGridOrbitals& phi, LocGridOrbitals& res, bool xlbomd);

    double evaluateEnergy(const LocGridOrbitals& phi);
};

#endif
