// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_MVP_DMStrategy_H
#define MGMOL_MVP_DMStrategy_H

#include "DMStrategy.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "MGmol.h"
#include "ProjectedMatricesInterface.h"
#include "Rho.h"

#include <iostream>
#include <mpi.h>
#include <vector>

template <class OrbitalsType, class MatrixType>
class MVP_DMStrategy : public DMStrategy<OrbitalsType>
{
private:
    ProjectedMatricesInterface* proj_matrices_;

    MPI_Comm comm_;
    std::ostream& os_;

    Ions& ions_;
    Rho<OrbitalsType>* rho_;
    Energy<OrbitalsType>* energy_;
    Electrostatic* electrostat_;
    Hamiltonian<OrbitalsType>* hamiltonian_;
    const std::vector<std::vector<int>>& global_indexes_;
    MGmol<OrbitalsType>* mgmol_strategy_;

    bool use_old_dm_;

public:
    MVP_DMStrategy(MPI_Comm comm, std::ostream& os, Ions& ions,
        Rho<OrbitalsType>* rho, Energy<OrbitalsType>* energy,
        Electrostatic* electrostat, Hamiltonian<OrbitalsType>*,
        MGmol<OrbitalsType>* mgmol_strategy,
        const std::vector<std::vector<int>>& overlappingGids,
        ProjectedMatricesInterface* proj_matrices, const bool use_old_dm);

    void initialize(OrbitalsType&) override {};
    int update(OrbitalsType& orbitals) override;

    // H is updated with MVP loop, so no need to compute it outside
    bool needH() const override { return false; }

    void stripDM() override;

    void dressDM() override;

    void reset() override { }
};

#endif
