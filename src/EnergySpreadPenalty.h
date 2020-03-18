// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef EnergySpreadPenalty_H
#define EnergySpreadPenalty_H

#include "SpreadPenaltyInterface.h"
#include "SpreadsAndCenters.h"

template <class T>
class EnergySpreadPenalty : public SpreadPenaltyInterface<T>
{
private:
    SpreadsAndCenters<T>* spreadf_;

    double spread2_target_;
    double alpha_;

    void computeAndAddResidualSpreadPenalty(
        const std::vector<float>& lagrangemult,
        const std::vector<float>& factors, const std::vector<Vector3D>& centers,
        const std::vector<int>& gids, T& orbitals, T& res);

public:
    EnergySpreadPenalty(SpreadsAndCenters<T>* spreadf,
        const double spread_target, const double alpha)
        : spreadf_(spreadf),
          spread2_target_(spread_target * spread_target),
          alpha_(alpha)
    {
        assert(spread2_target_ >= 0.);
        assert(alpha_ >= 0.);
    }

    // add penalty functional contribution to residual
    void addResidual(T& phi, T& res) override;

    double evaluateEnergy(const T& phi) override;
};

#endif
