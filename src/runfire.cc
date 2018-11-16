// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#include "Control.h"
#include "DFTsolver.h"
#include "FIRE.h"
#include "FIRE_IonicStepper.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "MPIdata.h"

#include <iomanip>
#include <iostream>
using namespace std;

void MGmol::runfire(LocGridOrbitals** orbitals, Ions& ions)
{
    printWithTimeStamp("Run FIRE algorithm...", cout);

    Control& ct = *(Control::instance());

    FIRE fire(orbitals, ions, *rho_, *constraints_, *lrs_, *currentMasks_,
        *electrostat_, ct.dt, *this);

    DFTsolver::resetItCount();

    fire.init(h5f_file_);

    delete h5f_file_;
    h5f_file_ = 0;

    // additional quench to compensate random start
    if (ct.restart_info < 3)
    {
        double eks = 0.;
        fire.quenchElectrons(ct.max_electronic_steps, eks);
    }
    else
    {
        DFTsolver::setItCountLarge();
    }

    // save computed vh for a fair energy "comparison" with vh computed
    // in close neigborhood

    // FIRE iterations
    for (int steps = 1; steps <= ct.num_MD_steps; steps++)
    {

        double eks;
        fire.quenchElectrons(ct.max_electronic_steps, eks);

        if (onpe0)
            os_ << setprecision(12) << fixed << "%%  " << steps
                << "  IONIC CONFIGURATION ENERGY = " << eks << endl;

        fire.computeForces();

        int flag_convF = fire.checkTolForces(ct.tol_forces);

        int conv = 0;
        if (flag_convF)
        {

            if (onpe0)
            {
                os_ << endl
                    << endl
                    << "FIRE: convergence in forces has been achieved. "
                       "stopping ..."
                    << endl;
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
                    << endl;
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
            if (onpe0) os_ << "Geometry optimization stopped" << endl;
            break;
        }

    } // end for steps

    // final dump
    if (ct.out_restart_info > 0)
    {
        fire.dumpRestart();
    }
}
