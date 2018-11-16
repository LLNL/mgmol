// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef ENERGY_H
#define ENERGY_H

#include "Grid.h"
#include "Timer.h"
#include "global.h"

#include <ostream>
#include <vector>

class Potentials;
class Ions;
class Electrostatic;
class LocGridOrbitals;
class ProjectedMatricesInterface;
class Rho;
class XConGrid;
class SpreadPenaltyInterface;

class Energy
{
    const pb::Grid& mygrid_;
    const Ions& ions_;
    const Potentials& pot_;
    const Electrostatic& es_;
    const Rho& rho_;
    const XConGrid& xc_;
    SpreadPenaltyInterface* spread_penalty_;

    std::vector<POTDTYPE> vofrho_;

    short nspin_;

    static Timer eval_te_tm_;

    double getEVrhoRho() const;

public:
    Energy(const pb::Grid&, const Ions&, const Potentials&,
        const Electrostatic&, const Rho&, const XConGrid&,
        SpreadPenaltyInterface*);

    static Timer eval_te_tm() { return eval_te_tm_; }

    double evaluateTotal(const double ts, ProjectedMatricesInterface*,
        const LocGridOrbitals& phi, const int, std::ostream&);

    double evaluateEnergyIonsInVext();

    void saveVofRho();
};

#endif
