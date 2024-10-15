// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef MGMOL_PoissonSolverFactory
#define MGMOL_PoissonSolverFactory

#include "Control.h"
#include "Hartree.h"
#include "Hartree_CG.h"
#include "Mesh.h"
#include "PBdiel.h"
#include "PBdiel_CG.h"
#include "ShiftedHartree.h"
#include "mputils.h"

#include "GridFunc.h"
#include "Laph2.h"
#include "Laph4.h"
#include "Laph4M.h"
#include "Laph4MP.h"
#include "Laph6.h"
#include "Laph8.h"
#include "ShiftedLaph4M.h"

class PoissonSolverFactory
{

public:
    /*!
     * return specific Poisson solver needed to solve Hartree problem
     */
    static Poisson* create(const pb::Grid& myGrid, PoissonFDtype lap_type,
        const short bc[3], const double screening_const)
    {
        Poisson* poisson_solver = nullptr;

        Control& ct = *(Control::instance());
        if (ct.MGPoissonSolver()) // use MG for Poisson Solver
        {
            if (screening_const > 0.)
            {
                switch (lap_type)
                {
                    case PoissonFDtype::h4M:
                        poisson_solver
                            = new ShiftedHartree<pb::ShiftedLaph4M<POTDTYPE>>(
                                myGrid, bc, screening_const);
                        break;
                    default:
                        (*MPIdata::sout)
                            << "Electrostatic, shifted, Undefined option: "
                            << static_cast<int>(lap_type) << std::endl;
                }
            }
            else
            {
                switch (lap_type)
                {
                    case PoissonFDtype::h4M:
                        poisson_solver
                            = new Hartree<pb::Laph4M<POTDTYPE>>(myGrid, bc);
                        break;
                    case PoissonFDtype::h2:
                        poisson_solver
                            = new Hartree<pb::Laph2<POTDTYPE>>(myGrid, bc);
                        break;
                    case PoissonFDtype::h4:
                        poisson_solver
                            = new Hartree<pb::Laph4<POTDTYPE>>(myGrid, bc);
                        break;
                    case PoissonFDtype::h6:
                        poisson_solver
                            = new Hartree<pb::Laph6<POTDTYPE>>(myGrid, bc);
                        break;
                    case PoissonFDtype::h8:
                        poisson_solver
                            = new Hartree<pb::Laph8<POTDTYPE>>(myGrid, bc);
                        break;
                    case PoissonFDtype::h4MP:
                        poisson_solver
                            = new Hartree<pb::Laph4MP<POTDTYPE>>(myGrid, bc);
                        break;
                    default:
                        (*MPIdata::sout)
                            << "Electrostatic, Undefined option: "
                            << static_cast<int>(lap_type) << std::endl;
                }
            }
        }
        else // use PCG for Poisson Solver
        {
            if (screening_const > 0.)
            {
                switch (lap_type)
                {
                    case PoissonFDtype::h4M:
                        poisson_solver
                            = new ShiftedHartree<pb::ShiftedLaph4M<POTDTYPE>>(
                                myGrid, bc, screening_const);
                        break;
                    default:
                        (*MPIdata::sout)
                            << "PCG Electrostatic, shifted, Undefined option: "
                            << static_cast<int>(lap_type) << std::endl;
                }
            }
            else
            {
                switch (lap_type)
                {
                    case PoissonFDtype::h4M:
                        poisson_solver
                            = new Hartree_CG<pb::Laph4M<POTDTYPE>>(myGrid, bc);
                        break;
                    case PoissonFDtype::h2:
                        poisson_solver
                            = new Hartree_CG<pb::Laph2<POTDTYPE>>(myGrid, bc);
                        break;
                    case PoissonFDtype::h4:
                        poisson_solver
                            = new Hartree_CG<pb::Laph4<POTDTYPE>>(myGrid, bc);
                        break;
                    case PoissonFDtype::h6:
                        poisson_solver
                            = new Hartree_CG<pb::Laph6<POTDTYPE>>(myGrid, bc);
                        break;
                    case PoissonFDtype::h8:
                        poisson_solver
                            = new Hartree_CG<pb::Laph8<POTDTYPE>>(myGrid, bc);
                        break;
                    case PoissonFDtype::h4MP:
                        poisson_solver
                            = new Hartree_CG<pb::Laph4MP<POTDTYPE>>(myGrid, bc);
                        break;
                    default:
                        (*MPIdata::sout)
                            << "PCG Electrostatic, Undefined option: "
                            << static_cast<int>(lap_type) << std::endl;
                }
            }
        }

        return poisson_solver;
    }

    static Poisson* createDiel(pb::Grid& pbGrid, PoissonFDtype lap_type,
        const short bc[3], const double e0, const double rho0,
        const double drho0)
    {
        Poisson* poisson_solver = nullptr;

        Control& ct = *(Control::instance());
        if (ct.MGPoissonSolver()) // use MG for Poisson Solver
        {
            switch (lap_type)
            {
                case PoissonFDtype::h4M:
                    poisson_solver = new PBdiel<pb::PBh4M<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                case PoissonFDtype::h2:
                    poisson_solver = new PBdiel<pb::PBh2<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                case PoissonFDtype::h4:
                    poisson_solver = new PBdiel<pb::PBh4<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                case PoissonFDtype::h6:
                    poisson_solver = new PBdiel<pb::PBh6<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                case PoissonFDtype::h8:
                    poisson_solver = new PBdiel<pb::PBh8<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                case PoissonFDtype::h4MP:
                    poisson_solver = new PBdiel<pb::PBh4MP<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                default:
                    (*MPIdata::sout)
                        << "Electrostatic, Undefined option" << std::endl;
            }
        }
        else // use PCG for Poisson Solver
        {
            switch (lap_type)
            {
                case PoissonFDtype::h4M:
                    poisson_solver = new PBdiel_CG<pb::PBh4M<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                case PoissonFDtype::h2:
                    poisson_solver = new PBdiel_CG<pb::PBh2<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                case PoissonFDtype::h4:
                    poisson_solver = new PBdiel_CG<pb::PBh4<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                case PoissonFDtype::h6:
                    poisson_solver = new PBdiel_CG<pb::PBh6<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                case PoissonFDtype::h8:
                    poisson_solver = new PBdiel_CG<pb::PBh8<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                case PoissonFDtype::h4MP:
                    poisson_solver = new PBdiel_CG<pb::PBh4MP<POTDTYPE>>(
                        pbGrid, bc, e0, rho0, drho0);
                    break;
                default:
                    (*MPIdata::sout)
                        << "Electrostatic, Undefined option" << std::endl;
            }
        }
        return poisson_solver;
    }
};

#endif
