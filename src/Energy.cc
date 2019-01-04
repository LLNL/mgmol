// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Energy.h"
#include "Control.h"
#include "Electrostatic.h"
#include "Ions.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"
#include "Potentials.h"
#include "ProjectedMatricesInterface.h"
#include "Rho.h"
#include "SpreadPenalty.h"
#include "XConGrid.h"

#include "Grid.h"
using namespace std;

#define RY2HA 0.5

template <class T>
Timer Energy<T>::eval_te_tm_("Energy::eval_te");

template <class T>
Energy<T>::Energy(const pb::Grid& mygrid, const Ions& ions, const Potentials& pot,
    const Electrostatic& es, const Rho<T>& rho,
    const XConGrid& xc,
    SpreadPenaltyInterface<T>* spread_penalty)
    : mygrid_(mygrid),
      ions_(ions),
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
    if (onpe0) (*MPIdata::sout) << "Energy<T>::saveVofRho()" << endl;
#endif
    pot_.getVofRho(vofrho_);
}

// get integral of rho*v[rho] in Hartree
template <class T>
double Energy<T>::getEVrhoRho() const
{
    double e = rho_.dotWithRho(&vofrho_[0]);

    // if(onpe0)(*MPIdata::sout)<<"get_Evrhorho="<<e<<endl;
    return RY2HA * mygrid_.vel() * e;
}

template <class T>
double Energy<T>::evaluateEnergyIonsInVext()
{
    double energy = 0.;

#ifdef HAVE_TRICUBIC
    if (!pot_.withVext()) return energy;

    //(*MPIdata::sout)<<"Energy<T>::evaluateEnergyIonsInVext()"<<endl;
    double position[3];
    vector<double> positions;
    positions.reserve(3 * ions_.local_ions().size());

    // loop over ions
    int nions                        = 0;
    vector<Ion*>::const_iterator ion = ions_.local_ions().begin();
    while (ion != ions_.local_ions().end())
    {
        (*ion)->getPosition(position);
        positions.push_back(position[0]);
        positions.push_back(position[1]);
        positions.push_back(position[2]);
        nions++;
        ion++;
    }

    vector<double> val(nions);
    pot_.getValVext(positions, val);

    // loop over ions again
    ion           = ions_.local_ions().begin();
    int ion_index = 0;
    while (ion != ions_.local_ions().end())
    {
        const double z = (*ion)->getZion();
        // int ion_index=(*ion)->index();
        assert(ion_index < (int)val.size());
        energy -= val[ion_index] * z;
        ion++;
        ion_index++;
    }

    //(*MPIdata::sout)<<"Energy<T>::evaluateEnergyIonsInVext(),
    //energy="<<energy<<endl;
    double tmp      = 0.;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&energy, &tmp, 1, MPI_SUM);
    energy = tmp;
#endif
    return energy;
}

template <class T>
double Energy<T>::evaluateTotal(const double ts, // in [Ha]
    ProjectedMatricesInterface* projmatrices, const T& phi,
    const int verbosity, ostream& os)
{
    eval_te_tm_.start();

    Control& ct = *(Control::instance());

    const double eself = ions_.energySelf();
    const double ediff = ions_.energyDiff(ct.bcPoisson);
    const double eipot = evaluateEnergyIonsInVext();

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
        os << "E epsilon=" << eepsilon << endl;
#endif
    }

    double spread_penalty_energy = 0.;
    if (spread_penalty_ != 0)
    {
        spread_penalty_energy = 0.5 * spread_penalty_->evaluateEnergy(phi);
        energy_sc += spread_penalty_energy;
    }

    if (verbosity > 1)
    {
        const vector<POTDTYPE>& vnuc(pot_.vnuc()); // vnuc is in [Ha]
        double evnuc = rho_.dotWithRho(&vnuc[0]);
        evnuc        = mygrid_.vel() * evnuc;
        if (onpe0)
        {
            os << setprecision(8) << fixed << endl;
            if (pot_.diel()) os << " E epsilon      [Ha] =" << eepsilon << endl;
            os << " EIGENVALUE SUM [Ha] = " << eigsum << endl;
            if (verbosity > 2)
            { // terms independent of electronic structure
                os << " Ions Eself     [Ha] = " << eself << endl;
                os << " Ions Ediff     [Ha] = " << ediff << endl;
            }
            os << " Ees+Eps+Eii    [Ha] = " << ees + evnuc - eself + ediff
               << endl;
            os << " E VH*RHO       [Ha] = " << Evh_rho << endl;
            os << " V*RHO          [Ha] = " << evrho << endl;
            os << " V_l*Rho        [Ha] = " << evnuc << endl;
            os << " ELECTROSTATIC  [Ha] = " << ees << endl;
            os << " XC             [Ha] = " << exc << endl;
            os << " -TS            [Ha] = " << -ts << endl;
            if (spread_penalty_ != 0)
                os << " Spread Penalty [Ha] = " << spread_penalty_energy
                   << endl;
            os << " Ions Ext. Pot. [Ha] = " << eipot << endl;

            os << " SC ENERGY      [Ha] = " << energy_sc << endl;
        }
    }

    eval_te_tm_.stop();

    return energy_sc;
}

template class Energy<LocGridOrbitals>;
template class Energy<ExtendedGridOrbitals>;

