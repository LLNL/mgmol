// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_SPREADPENALTYVOLUME_H
#define MGMOL_SPREADPENALTYVOLUME_H

#include "SpreadPenaltyInterface.h"
#include "SpreadsAndCenters.h"

template <class T>
class SpreadPenaltyVolume : public SpreadPenaltyInterface<T>
{
private:
    SpreadsAndCenters<T>* spreadf_;

    double spread_target_;
    double alpha_;
    double dampingFactor_;

    double volume_target_;

    void computeAndAddResidualSpreadPenalty(
        const std::vector<float>& lagrangemult, const float factor,
        const std::vector<Vector3D>& centers, const std::vector<int>& gids,
        T& orbitals, T& res);

public:
    SpreadPenaltyVolume(SpreadsAndCenters<T>* spreadf,
        const double spread_target, const double alpha,
        const double dampingFactor)
        : spreadf_(spreadf),
          spread_target_(spread_target),
          alpha_(alpha),
          dampingFactor_(dampingFactor)
    {
        assert(spread_target_ >= 0.);
        assert(alpha_ >= 0.);
        assert(dampingFactor_ > 0.);

        volume_target_ = (spread_target_ * spread_target_ * spread_target_);
    }

    // add penalty functional contribution to residual
    void addResidual(T& phi, T& res) override;
    double evaluateEnergy(const T& phi) override;
};

#endif
