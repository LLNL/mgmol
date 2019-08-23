// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DMSTRATEGYFACTORY_H
#define MGMOL_DMSTRATEGYFACTORY_H

#include "Control.h"
#include "DistMatrix.h"
#include "DistMatrixWithSparseComponent.h"
#include "EigenDMStrategy.h"
#include "FullyOccupiedNonOrthoDMStrategy.h"
#include "HamiltonianMVP_DMStrategy.h"
#include "MGmol.h"
#include "MVP_DMStrategy.h"
#include "NonOrthoDMStrategy.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesSparse.h"

template <class T>
class DMStrategyFactory
{
public:
    static DMStrategy* create(MPI_Comm comm, std::ostream& os, Ions& ions,
        Rho<T>* rho, Energy<T>* energy, Electrostatic* electrostat,
        MGmol<T>* mgmol_strategy, ProjectedMatricesInterface* proj_matrices,
        T* orbitals)
    {
        Control& ct = *(Control::instance());

        DMStrategy* dm_strategy;
        if (ct.DM_solver() == DMNonLinearSolverType::MVP)
        {
            dm_strategy = new MVP_DMStrategy<T>(comm, os, ions, rho, energy,
                electrostat, mgmol_strategy, orbitals, proj_matrices,
                ct.use_old_dm());
        }
        else if (ct.DM_solver() == DMNonLinearSolverType::HMVP)
        {
            dm_strategy = createHamiltonianMVP_DMStrategy(comm, os, ions, rho,
                energy, electrostat, mgmol_strategy, proj_matrices, orbitals,
                ct.short_sighted);
        }
        else
        {
            if (ct.fullyOccupied())
            {
                dm_strategy
                    = new FullyOccupiedNonOrthoDMStrategy(proj_matrices);
            }
            else
            {
                if (ct.getOrbitalsType() == OrbitalsType::Eigenfunctions)
                {
                    dm_strategy
                        = new EigenDMStrategy<T>(orbitals, proj_matrices);
                }
                else
                {
                    if (ct.getOrbitalsType() == OrbitalsType::Nonorthogonal)
                    {
                        dm_strategy = new NonOrthoDMStrategy<T>(
                            orbitals, proj_matrices, ct.dm_mix);
                    }
                }
            }
        }

        return dm_strategy;
    }

private:
    static DMStrategy* createHamiltonianMVP_DMStrategy(MPI_Comm comm,
        std::ostream& os, Ions& ions, Rho<T>* rho, Energy<T>* energy,
        Electrostatic* electrostat, MGmol<T>* mgmol_strategy,
        ProjectedMatricesInterface* proj_matrices, T*, const bool);
};

#endif
