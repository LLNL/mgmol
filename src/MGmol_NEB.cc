// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DFTsolver.h"
#include "FIRE.h"
#include "LBFGS.h"
#include "MGmol.h"
#include "ProjectedMatricesInterface.h"

using namespace std;

template <class T>
void MGmol<T>::sebprintForces()
{
    ions_->printForces(os_);
}
template <class T>
void MGmol<T>::sebprintPositions()
{
    ions_->printPositions(os_);
}

template <class T>
void MGmol<T>::geomOptimSetup()
{
    Control& ct = *(Control::instance());

    switch (ct.AtomsDynamic())
    {
        case AtomsDynamicType::LBFGS:
            geom_optimizer_ = new LBFGS<T>(&current_orbitals_, *ions_, *rho_,
                *constraints_, *lrs_, local_cluster_, *currentMasks_,
                *corrMasks_, *electrostat_, ct.dt, *this);
            break;

        case AtomsDynamicType::FIRE:
            geom_optimizer_
                = new FIRE<T>(&current_orbitals_, *ions_, *rho_, *constraints_,
                    *lrs_, *currentMasks_, *electrostat_, ct.dt, *this);
            break;

        default:
            (*MPIdata::serr) << "geomOptimSetup(): option "
                             << " is an invalid method" << endl;
            return;
    }
    DFTsolver<T>::resetItCount();

    geom_optimizer_->init(h5f_file_);

    // additional quench to compensate random start
    if (ct.restart_info < 3)
    {
        double eks = 0.;
        geom_optimizer_->quenchElectrons(ct.max_electronic_steps, eks);
    }
    else
    {
        DFTsolver<T>::setItCountLarge();
    }
}

template <class T>
void MGmol<T>::geomOptimQuench()
{
    Control& ct = *(Control::instance());

    double eks = 0.;
    geom_optimizer_->quenchElectrons(ct.max_electronic_steps, eks);

    // Get the total energy
    double ts     = 0.5 * proj_matrices_->computeEntropy(); // in [Ha]
    total_energy_ = energy_->evaluateTotal(
        ts, proj_matrices_, *current_orbitals_, 2, os_);
}

template <class T>
void MGmol<T>::geomOptimComputeForces()
{
    geom_optimizer_->computeForces();
}

template <class T>
void MGmol<T>::geomOptimSetForces(const vector<vector<double>>& f)
{
    geom_optimizer_->setForces(f);
}

template <class T>
void MGmol<T>::geomOptimDumpRestart()
{
    geom_optimizer_->dumpRestart();
}

template <class T>
int MGmol<T>::geomOptimRun1Step()
{
    int conv = geom_optimizer_->run1step();
    geom_optimizer_->updatePotAndMasks();
    // Write down positions and displacements
    ions_->printPositions(os_);
    return conv;
}

template <class T>
short MGmol<T>::geomOptimCheckTolForces(const double tol_force)
{
    return geom_optimizer_->checkTolForces(tol_force);
}
