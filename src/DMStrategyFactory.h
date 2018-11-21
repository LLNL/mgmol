// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef DMSTRATEGYFACTORY_H
#define DMSTRATEGYFACTORY_H

#include "Control.h"
#include "DistMatrix.h"
#include "DistMatrixWithSparseComponent.h"
#include "EigenDMStrategy.h"
#include "FullyOccupiedNonOrthoDMStrategy.h"
#include "HamiltonianMVP_DMStrategy.h"
#include "MVP_DMStrategy.h"
#include "NonOrthoDMStrategy.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesSparse.h"

class DMStrategyFactory
{
public:
    static DMStrategy* create(MPI_Comm comm, std::ostream& os, Ions& ions,
        Rho* rho, Energy* energy, Electrostatic* electrostat,
        MGmol* mgmol_strategy, ProjectedMatricesInterface* proj_matrices,
        LocGridOrbitals* orbitals)
    {
        Control& ct = *(Control::instance());

        DMStrategy* dm_strategy;
        if (ct.DM_solver() == DMNonLinearSolverType::MVP)
        {
            dm_strategy
                = new MVP_DMStrategy(comm, os, ions, rho, energy, electrostat,
                    mgmol_strategy, orbitals, proj_matrices, ct.use_old_dm());
        }
        else if (ct.DM_solver() == DMNonLinearSolverType::HMVP)
        {
            if (ct.short_sighted)
            {
                dm_strategy = new HamiltonianMVP_DMStrategy<
                    VariableSizeMatrix<sparserow>,
                    VariableSizeMatrix<sparserow>, ProjectedMatricesSparse>(
                    comm, os, ions, rho, energy, electrostat, mgmol_strategy,
                    orbitals);
            }
            else
            {
                dm_strategy = new HamiltonianMVP_DMStrategy<
                    dist_matrix::DistMatrix<DISTMATDTYPE>,
                    dist_matrix::DistMatrixWithSparseComponent<DISTMATDTYPE>,
                    ProjectedMatrices>(comm, os, ions, rho, energy, electrostat,
                    mgmol_strategy, orbitals);
            }
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
                    dm_strategy = new EigenDMStrategy(orbitals, proj_matrices);
                }
                else
                {
                    if (ct.getOrbitalsType() == OrbitalsType::Nonorthogonal)
                    {
                        dm_strategy = new NonOrthoDMStrategy(
                            orbitals, proj_matrices, ct.dm_mix);
                    }
                }
            }
        }

        return dm_strategy;
    }
};

#endif
