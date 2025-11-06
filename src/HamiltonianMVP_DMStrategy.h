// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_HamiltonianMVP_DMStrategy_H
#define MGMOL_HamiltonianMVP_DMStrategy_H

#include "DMStrategy.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "HamiltonianMVPSolver.h"
#include "Ions.h"
#include "MGmol.h"
#include "Rho.h"

template <class T>
class MGmol;

template <class MatrixType, class ProjMatrixType, class OrbitalsType>
class HamiltonianMVP_DMStrategy : public DMStrategy<OrbitalsType>
{
private:
    MPI_Comm comm_;
    std::ostream& os_;

    Ions& ions_;
    Rho<OrbitalsType>* rho_;
    Energy<OrbitalsType>* energy_;
    Electrostatic* electrostat_;
    Hamiltonian<OrbitalsType>* hamiltonian_;
    const std::vector<std::vector<int>>& global_indexes_;
    MGmol<OrbitalsType>* mgmol_strategy_;

    HamiltonianMVPSolver<MatrixType, ProjMatrixType, OrbitalsType>* solver_;

public:
    HamiltonianMVP_DMStrategy(MPI_Comm comm, std::ostream& os, Ions& ions,
        Rho<OrbitalsType>* rho, Energy<OrbitalsType>* energy,
        Electrostatic* electrostat, Hamiltonian<OrbitalsType>* hamiltonian,
        MGmol<OrbitalsType>* mgmol_strategy, OrbitalsType* orbitals);

    ~HamiltonianMVP_DMStrategy() override;

    void initialize(OrbitalsType& orbitals) override;
    int update(OrbitalsType& orbitals) override;

    // H is updated with HamiltonianMVP loop, so no need to compute it outside
    bool needH() const override { return false; }

    void stripDM() override;

    void dressDM() override;

    void reset() override;
};

#endif
