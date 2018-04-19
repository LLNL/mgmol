// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef SpreadPenaltyVolume_H
#define SpreadPenaltyVolume_H

class LocGridOrbitals;
#include "SpreadsAndCenters.h"
#include "SpreadPenaltyInterface.h"

class SpreadPenaltyVolume: public SpreadPenaltyInterface
{
private:

    SpreadsAndCenters* spreadf_;

    double spread_target_;
    double alpha_;
    double dampingFactor_;
    
    double volume_target_;
    
    void computeAndAddResidualSpreadPenalty(const vector<float>& lagrangemult,
                                            const float factor,
                                            const vector<Vector3D>& centers,
                                            const vector<int>& gids,
                                            LocGridOrbitals& orbitals,
                                            LocGridOrbitals& res);
public:

    SpreadPenaltyVolume(SpreadsAndCenters* spreadf,
                  const double spread_target,
                  const double alpha,
                  const double dampingFactor):
        spreadf_(spreadf),
        spread_target_(spread_target),
        alpha_(alpha),
        dampingFactor_(dampingFactor)
    {
        assert( spread_target_>=0. );
        assert( alpha_>=0. );
        assert( dampingFactor_>0. );
        
        volume_target_=(spread_target_*spread_target_*spread_target_);
    }

    //add penalty functional contribution to residual
    void addResidual(LocGridOrbitals& phi,
                     LocGridOrbitals& res);
    double evaluateEnergy(const LocGridOrbitals& phi);
};

#endif
