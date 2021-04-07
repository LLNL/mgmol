// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DielectricControl.h"

#include "Control.h"

DielectricControl::DielectricControl(Potentials& pot,
    Electrostatic* electrostat, std::ostream& os, const bool onpe0)
    : pot_(pot),
      electrostat_(electrostat),
      os_(os),
      onpe0_(onpe0),
      pbset_(false)
{
}

void DielectricControl::activate(const int it_scf, const double deig2)
{
    bool isON = pot_.diel();
    if (!isON) return; // continuum solvent is OFF

    if (pbset_) return; // continuum solvent already set

    Control& ct(*(Control::instance()));
    const int diel_delay = (ct.restart_info < 3) ? 10 : 1;
    if (it_scf < diel_delay)
    {
        isON = false;
    }

    // turn ON continuum solvent
    if (isON && (ct.restart_info >= 3 || deig2 < 2.e-2 * ct.numst))
    {
        electrostat_->setupPB(ct.e0_, ct.rho0_, ct.drho0_, pot_);
        electrostat_->setup(ct.vh_its);
        pbset_ = true;
    }

    if (pot_.diel() && !pbset_)
        if (onpe0_ && ct.verbose > 1)
            os_ << " Solvation turned off for this step" << std::endl;
}
