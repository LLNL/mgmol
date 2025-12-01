// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Energy.h"
#include "Control.h"
#include "Electrostatic.h"
#include "ExtendedGridOrbitals.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"
#include "Potentials.h"
#include "ProjectedMatricesInterface.h"
#include "Rho.h"
#include "SpreadPenalty.h"
#include "XConGrid.h"

#include "Grid.h"

#define RY2HA 0.5

template <class T>
Energy<T>::Energy(const pb::Grid& mygrid, const Potentials& pot,
    const Electrostatic& es, const Rho<T>& rho, const XConGrid& xc,
    SpreadPenaltyInterface<T>* spread_penalty)
    : mygrid_(mygrid),
      pot_(pot),
      es_(es),
      rho_(rho),
      xc_(xc),
      spread_penalty_(spread_penalty)
{
    nspin_ = rho.rho_.size();
}

template <class T>
void Energy<T>::saveVofRho()
{
#ifdef PRINT_OPERATIONS
    if (onpe0) (*MPIdata::sout) << "Energy<T>::saveVofRho()" << std::endl;
#endif
    pot_.getVofRho(vofrho_);
}

// get integral of rho*v[rho] in Hartree
template <class T>
double Energy<T>::getEVrhoRho() const
{
    double e = rho_.dotWithRho(&vofrho_[0]);

    // if(onpe0)(*MPIdata::sout)<<"get_Evrhorho="<<e<<std::endl;
    return RY2HA * mygrid_.vel() * e;
}

template <class T>
double Energy<T>::evaluateEnergyIonsInVext(Ions& ions)
{
#ifdef HAVE_TRICUBIC
    double energy = 0.;

    if (!pot_.withVext()) return energy;

    //(*MPIdata::sout)<<"Energy<T>::evaluateEnergyIonsInVext()"<<std::endl;
    double position[3];
    std::vector<double> positions;
    positions.reserve(3 * ions.local_ions().size());

    // loop over ions
    int nions                             = 0;
    std::vector<Ion*>::const_iterator ion = ions.local_ions().begin();
    while (ion != ions.local_ions().end())
    {
        (*ion)->getPosition(position);
        positions.push_back(position[0]);
        positions.push_back(position[1]);
        positions.push_back(position[2]);
        nions++;
        ion++;
    }

    std::vector<double> val(nions);
    pot_.getValVext(positions, val);

    double energy = 0.;

    // loop over ions again
    ion           = ions.local_ions().begin();
    int ion_index = 0;
    while (ion != ions.local_ions().end())
    {
        const double z = (*ion)->getZion();
        // int ion_index=(*ion)->index();
        assert(ion_index < (int)val.size());
        energy -= val[ion_index] * z;
        ion++;
        ion_index++;
    }

    //(*MPIdata::sout)<<"Energy<T>::evaluateEnergyIonsInVext(),
    // energy="<<energy<<std::endl;
    double tmp      = 0.;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&energy, &tmp, 1, MPI_SUM);
    energy = tmp;

    return energy;
#else
    (void)ions;

    return 0.;
#endif
}

template <class T>
double Energy<T>::evaluateTotal(const double ts, // in [Ha]
    ProjectedMatricesInterface* projmatrices, Ions& ions, const T& phi,
    const int verbosity, std::ostream& os)
{
    eval_te_tm_.start();

    Control& ct = *(Control::instance());

    const double eself = ions.energySelf();
    const double ediff = ions.energyDiff(ct.bcPoisson);
    const double eipot = evaluateEnergyIonsInVext(ions);

    const double eigsum = 0.5 * projmatrices->getExpectationH();

    // electrostatic energy
    const double Evh_rho  = es_.evhRho();
    const double Evh_rhoc = es_.evhRhoc();
    const double ees      = 0.5 * (Evh_rho - Evh_rhoc);

    const double exc = xc_.getExc();

    const double evrho = getEVrhoRho();

    double energy_sc = eigsum - evrho + exc + ees - eself + ediff + eipot - ts;
    // Correct energy:
    // substract contribution Vespilon in sum over eigenvalues
    const double eepsilon = es_.eEpsilon();
    if (pot_.diel())
    {
        energy_sc -= eepsilon;
#ifdef DEBUG
        os << "E epsilon=" << eepsilon << std::endl;
#endif
    }

    double spread_penalty_energy = 0.;
    if (spread_penalty_ != nullptr)
    {
        spread_penalty_energy = 0.5 * spread_penalty_->evaluateEnergy(phi);
        energy_sc += spread_penalty_energy;
    }

    if (verbosity > 1)
    {
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());

        const std::vector<POTDTYPE>& vnuc(pot_.vnuc()); // vnuc is in [Ha]
        double evnuc = rho_.dotWithRho(&vnuc[0]);
        evnuc        = mygrid_.vel() * evnuc;
        if (mmpi.PE0())
        {
            os << std::setprecision(8) << std::fixed << std::endl;
            if (pot_.diel())
                os << " E epsilon      [Ha] =" << eepsilon << std::endl;
            os << " EIGENVALUE SUM [Ha] = " << eigsum << std::endl;
            if (verbosity > 2)
            { // terms independent of electronic structure
                os << " Ions Eself     [Ha] = " << eself << std::endl;
                os << " Ions Ediff     [Ha] = " << ediff << std::endl;
            }
            os << " Ees+Eps+Eii    [Ha] = " << ees + evnuc - eself + ediff
               << std::endl;
            os << " E VH*RHO       [Ha] = " << Evh_rho << std::endl;
            os << " V*RHO          [Ha] = " << evrho << std::endl;
            os << " V_l*Rho        [Ha] = " << evnuc << std::endl;
            os << " ELECTROSTATIC  [Ha] = " << ees << std::endl;
            os << " XC             [Ha] = " << exc << std::endl;
            os << " -TS            [Ha] = " << -ts << std::endl;
            if (spread_penalty_ != nullptr)
                os << " Spread Penalty [Ha] = " << spread_penalty_energy
                   << std::endl;
            os << " Ions Ext. Pot. [Ha] = " << eipot << std::endl;

            os << " SC ENERGY      [Ha] = " << energy_sc << std::endl;
        }
    }

    eval_te_tm_.stop();

    return energy_sc;
}

template class Energy<LocGridOrbitals<ORBDTYPE>>;
template class Energy<ExtendedGridOrbitals<ORBDTYPE>>;
