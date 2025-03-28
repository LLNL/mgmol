// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ENERGY_H
#define MGMOL_ENERGY_H

#include "Grid.h"
#include "Ions.h"
#include "Rho.h"
#include "SpreadPenaltyInterface.h"
#include "Timer.h"
#include "global.h"

#include <ostream>
#include <vector>

class Potentials;
class Electrostatic;
class ProjectedMatricesInterface;
class XConGrid;

template <class T>
class Energy
{
    const pb::Grid& mygrid_;
    const Potentials& pot_;
    const Electrostatic& es_;
    const Rho<T>& rho_;
    const XConGrid& xc_;
    SpreadPenaltyInterface<T>* spread_penalty_;

    std::vector<POTDTYPE> vofrho_;

    short nspin_;

    static Timer eval_te_tm_;

    double getEVrhoRho() const;

public:
    Energy(const pb::Grid&, const Potentials&, const Electrostatic&,
        const Rho<T>&, const XConGrid&, SpreadPenaltyInterface<T>*);

    static Timer eval_te_tm() { return eval_te_tm_; }

    double evaluateTotal(const double ts, ProjectedMatricesInterface*, Ions&,
        const T& phi, const int, std::ostream&);

    double evaluateEnergyIonsInVext(Ions&);

    void saveVofRho();
};
// Instantiate static variable here to avoid clang warnings
template <class T>
Timer Energy<T>::eval_te_tm_("Energy::eval_te");
#endif
