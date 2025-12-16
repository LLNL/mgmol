// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DMSTRATEGYFACTORY_H
#define MGMOL_DMSTRATEGYFACTORY_H

#include "Control.h"
#include "EigenDMStrategy.h"
#include "FullyOccupiedNonOrthoDMStrategy.h"
#include "HamiltonianMVP_DMStrategy.h"
#include "MGmol.h"
#include "MVP_DMStrategy.h"
#include "NonOrthoDMStrategy.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesSparse.h"

template <class OrbitalsType, class MatrixType>
class DMStrategyFactory
{
public:
    static DMStrategy<OrbitalsType>* create(MPI_Comm comm, std::ostream& os,
        Ions& ions, Rho<OrbitalsType>* rho, Energy<OrbitalsType>* energy,
        Electrostatic* electrostat, Hamiltonian<OrbitalsType>* hamiltonian,
        MGmol<OrbitalsType>* mgmol_strategy,
        ProjectedMatricesInterface* proj_matrices, OrbitalsType* orbitals)
    {
        Control& ct     = *(Control::instance());
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());

        DMStrategy<OrbitalsType>* dm_strategy = nullptr;
        if (ct.DM_solver() == DMNonLinearSolverType::MVP)
        {
            dm_strategy = new MVP_DMStrategy<OrbitalsType, MatrixType>(comm, os,
                ions, rho, energy, electrostat, hamiltonian, mgmol_strategy,
                orbitals->getOverlappingGids(), proj_matrices, ct.use_old_dm());
        }
        else if (ct.DM_solver() == DMNonLinearSolverType::HMVP)
        {
            dm_strategy = createHamiltonianMVP_DMStrategy(comm, os, ions, rho,
                energy, electrostat, hamiltonian, mgmol_strategy, proj_matrices,
                orbitals, ct.short_sighted);
        }
        else
        {
            if (ct.fullyOccupied())
            {
                if (mmpi.instancePE0())
                    std::cout << "Fully occupied strategy" << std::endl;
                dm_strategy = new FullyOccupiedNonOrthoDMStrategy<OrbitalsType>(
                    proj_matrices);
            }
            else
            {
                if (ct.getOrthoType() == OrthoType::Eigenfunctions)
                {
                    if (mmpi.instancePE0())
                        std::cout << "EigenDMStrategy..." << std::endl;
                    dm_strategy
                        = new EigenDMStrategy<OrbitalsType>(proj_matrices);
                }
                else
                {
                    if (ct.getOrthoType() == OrthoType::Nonorthogonal)
                    {
                        if (mmpi.instancePE0())
                            std::cout << "NonOrthoDMStrategy..." << std::endl;
                        dm_strategy = new NonOrthoDMStrategy<OrbitalsType>(
                            proj_matrices, ct.dm_mix);
                    }
                }
            }
        }

        assert(dm_strategy != nullptr);
        return dm_strategy;
    }

private:
    static DMStrategy<OrbitalsType>* createHamiltonianMVP_DMStrategy(
        MPI_Comm comm, std::ostream& os, Ions& ions, Rho<OrbitalsType>* rho,
        Energy<OrbitalsType>* energy, Electrostatic* electrostat,
        Hamiltonian<OrbitalsType>* hamiltonian,
        MGmol<OrbitalsType>* mgmol_strategy,
        ProjectedMatricesInterface* proj_matrices, OrbitalsType* orbitals,
        const bool);
};

#endif
