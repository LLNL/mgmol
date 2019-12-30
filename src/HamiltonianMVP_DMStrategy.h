// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_HamiltonianMVP_DMStrategy_H
#define MGMOL_HamiltonianMVP_DMStrategy_H

#include "DMStrategy.h"
#include "Energy.h"
#include "HamiltonianMVPSolver.h"
#include "MGmol.h"
#include "Rho.h"

class Ions;
class Electrostatic;
template <class T>
class MGmol;

template <class T1, class T3, class T4>
class HamiltonianMVP_DMStrategy : public DMStrategy
{
private:
    T4* orbitals_;

    MPI_Comm comm_;
    std::ostream& os_;

    Ions& ions_;
    Rho<T4>* rho_;
    Energy<T4>* energy_;
    Electrostatic* electrostat_;
    const std::vector<std::vector<int>>& global_indexes_;
    MGmol<T4>* mgmol_strategy_;

    HamiltonianMVPSolver<T1, T3, T4>* solver_;

public:
    HamiltonianMVP_DMStrategy(MPI_Comm comm, std::ostream& os, Ions& ions,
        Rho<T4>* rho, Energy<T4>* energy, Electrostatic* electrostat,
        MGmol<T4>* mgmol_strategy, T4* orbitals);

    ~HamiltonianMVP_DMStrategy() override;

    void initialize() override;
    int update() override;

    // H is updated with HamiltonianMVP loop, so no need to compute it outside
    bool needH() const override { return false; }

    void stripDM() override;

    void dressDM() override;

    void reset() override;
};

#endif
