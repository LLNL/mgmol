// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Control.h"
#include "DFTsolver.h"
#include "ExtendedGridOrbitals.h"
#include "FIRE.h"
#include "FIRE_IonicStepper.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "MPIdata.h"

#include <iomanip>
#include <iostream>

template <class OrbitalsType>
void MGmol<OrbitalsType>::runfire(OrbitalsType** orbitals, Ions& ions)
{
    printWithTimeStamp("Run FIRE algorithm...", std::cout);

    Control& ct = *(Control::instance());

    FIRE<OrbitalsType> fire(orbitals, ions, *rho_, *constraints_, lrs_,
        *currentMasks_, *electrostat_, ct.dt, *this);

    DFTsolver<OrbitalsType>::resetItCount();

    fire.init(h5f_file_);

    delete h5f_file_;
    h5f_file_ = nullptr;

    // additional quench to compensate random start
    if (ct.restart_info < 3)
    {
        double eks = 0.;
        fire.quenchElectrons(ct.max_electronic_steps, eks);
    }
    else
    {
        DFTsolver<OrbitalsType>::setItCountLarge();
    }

    // save computed vh for a fair energy "comparison" with vh computed
    // in close neigborhood

    // FIRE iterations
    for (int steps = 1; steps <= ct.num_MD_steps; steps++)
    {

        double eks;
        fire.quenchElectrons(ct.max_electronic_steps, eks);

        if (onpe0)
            os_ << std::setprecision(12) << std::fixed << "%%  " << steps
                << "  IONIC CONFIGURATION ENERGY = " << eks << std::endl;

        fire.computeForces();

        int flag_convF = fire.checkTolForces(ct.tol_forces);

        int conv = 0;
        if (flag_convF)
        {

            if (onpe0)
            {
                os_ << std::endl
                    << std::endl
                    << "FIRE: convergence in forces has been achieved. "
                       "stopping ..."
                    << std::endl;
            }
            conv = 1;
        }
        else
        {

            // tentative step for atomic positions
            conv = fire.run1step();

            if (onpe0)
            {
                os_ << "FIRE: update atomic configuration dependent stuff..."
                    << std::endl;
            }
            // update stuff that depends on atomic positions
            fire.updatePotAndMasks();

            if (ct.checkpoint && ct.out_restart_file != "0")
                if (ct.out_restart_info > 0)
                    if ((steps % ct.checkpoint) == 0 && steps < ct.num_MD_steps)
                    {
                        fire.dumpRestart();
                    }
        }

        // Write down positions and displacements
        if (onpe0) ions.printPositions(os_);

        if (conv != 0)
        {
            if (onpe0) os_ << "Geometry optimization stopped" << std::endl;
            break;
        }

    } // end for steps

    // final dump
    if (ct.out_restart_info > 0)
    {
        fire.dumpRestart();
    }
}

template void MGmol<LocGridOrbitals>::runfire(
    LocGridOrbitals** orbitals, Ions& ions);
template void MGmol<ExtendedGridOrbitals>::runfire(
    ExtendedGridOrbitals** orbitals, Ions& ions);
