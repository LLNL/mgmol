// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ConstraintSet.h"
#include "Control.h"
#include "DFTsolver.h"
#include "Electrostatic.h"
#include "Energy.h"
#include "Ions.h"
#include "LBFGS.h"
#include "LocalizationRegions.h"
#include "MGmol.h"
#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "MasksSet.h"
#include "Mesh.h"
#include "Potentials.h"
#include "ProjectedMatrices.h"
#include "Rho.h"

#include <iomanip>
#include <iostream>
#include <vector>
using namespace std;

template <class OrbitalsType>
void MGmol<OrbitalsType>::lbfgsrlx(OrbitalsType** orbitals, Ions& ions)
{
    Control& ct = *(Control::instance());

    LBFGS<OrbitalsType> lbfgs(orbitals, ions, *rho_, *constraints_, lrs_,
        local_cluster_.get(), *currentMasks_, *corrMasks_, *electrostat_, ct.dt,
        *this);

    DFTsolver<OrbitalsType>::resetItCount();

    lbfgs.init(h5f_file_.get());

    h5f_file_.reset();

    // additional quench to compensate random start
    if (ct.restart_info < 3)
    {
        double eks = 0.;
        lbfgs.quenchElectrons(ct.max_electronic_steps, eks);
    }
    else
    {
        DFTsolver<OrbitalsType>::setItCountLarge();
    }

    // save computed vh for a fair energy "comparison" with vh computed
    // in close neigborhood

    lbfgs.setupConstraints();

    // LBFGS iterations
    for (int steps = 1; steps <= ct.num_MD_steps; steps++)
    {
        if (steps > 3 || ct.restart_info > 2) lbfgs.setQuenchTol();

        if (steps > 1) lbfgs.updateRefs();

        double eks;
        lbfgs.quenchElectrons(ct.max_electronic_steps, eks);

        if (onpe0)
            os_ << setprecision(12) << fixed << "%%  " << steps
                << "  IONIC CONFIGURATION ENERGY = " << eks << endl;

        lbfgs.computeForces();

        // check for convergence
        int flag_convF = lbfgs.checkTolForces(ct.tol_forces);

        int conv = 0;
        if (flag_convF)
        {
            if (onpe0)
            {
                os_ << endl
                    << endl
                    << "LBFGS: convergence in forces has been achieved. "
                    << endl;
            }
        }
        else
        {
            // tentative step for atomic positions
            conv = lbfgs.run1step();

            if (onpe0)
            {
                os_ << "LBFGS: update atomic configuration dependent stuff..."
                    << endl;
            }
            // update stuff that depends on atomic positions
            lbfgs.updatePotAndMasks();

            if (ct.checkpoint && ct.out_restart_file != "0")
                if (ct.out_restart_info > 0)
                    if ((steps % ct.checkpoint) == 0 && steps < ct.num_MD_steps)
                    {
                        lbfgs.dumpRestart();
                    }
        }

        // Write down positions and displacements
        ions.printPositions(os_);

        // early termination if convergence achieved
        if (flag_convF)
        {
            if (onpe0) os_ << "LBFGS: stopping ..." << endl;
            break;
        }

        // reset iterations if last step failed
        if (conv != 0)
        {
            if (onpe0) os_ << "LBFGS Geometry optimization stopped" << endl;
            if (onpe0) os_ << "LBFGS Geometry optimization reset" << endl;
            lbfgs.reset(ct.dt);
        }

    } // end for steps

    // final dump
    if (ct.out_restart_info > 0)
    {
        lbfgs.dumpRestart();
    }
}

template void MGmol<LocGridOrbitals<ORBDTYPE>>::lbfgsrlx(
    LocGridOrbitals<ORBDTYPE>** orbitals, Ions& ions);
template void MGmol<ExtendedGridOrbitals<ORBDTYPE>>::lbfgsrlx(
    ExtendedGridOrbitals<ORBDTYPE>** orbitals, Ions& ions);
