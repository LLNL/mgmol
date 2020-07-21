// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#ifndef GlobalLBFGS_H
#define GlobalLBFGS_H

#include "ConstraintSet.h"
#include "Energy.h"
#include "GridFunc.h"
#include "Ions.h"
#include "LBFGS_IonicStepper.h"
#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "Rho.h"

#include <vector>

class MasksSet;
class Energy;
class Electrostatic;
class MGmol;
class KBPsiMatrixInterface;

class GlobalLBFGS
{
private:
    LocGridOrbitals** orbitals_;
    Ions& ions_;
    Rho& rho_;
    ConstraintSet& constraints_;
    KBPsiMatrixInterface* g_kbpsi_;
    LBFGS_IonicStepper* stepper_;
    Energy& energy_;
    std::shared_ptr<LocalizationRegions> lrs_;
    MasksSet& masks_;
    Electrostatic& electrostat_;
    LocalizationRegions ref_lrs_;
    MasksSet* ref_masks_;

    LocGridOrbitals* ref_orbitals_;
    pb::GridFunc<POTDTYPE>* vh_init_;

    std::vector<double> tau0_; // tau0[3*ia+j]
    std::vector<double> taup_; // taup[3*ia+j]
    std::vector<double> fion_; // fion[3*ia+j]
    std::vector<double> pmass_; // pmass[ia]
    std::vector<short> atmove_; // atmove[ia]
    std::vector<string> ions_names_;

    std::vector<double> gtau0_; // tau0[3*ia+j]
    std::vector<double> gtaup_; // taup[3*ia+j]
    std::vector<double> gfion_; // fion[3*ia+j]
    std::vector<short> gatmove_; // atmove[ia]

    double etot_i_[3];

    MGmol& mgmol_strategy_;
    int local_image_;

public:
    GlobalLBFGS(LocGridOrbitals** orbitals, Ions& ions, Rho& rho,
        ConstraintSet& constraints, KBPsiMatrixInterface* g_kbpsi,
        Energy& energy, std::shared_ptr<LocalizationRegions> lrs,
        MasksSet& masks, Electrostatic& electrostat, const double dt, MGmol&,
        const int local_image, const int nimages);

    ~GlobalLBFGS();

    void init(HDFrestart* h5f_file);
    int run1step();
    void computeForces();
    void printForces(ostream& os, const int root = 0);
    bool checkTolForces(const double tol);
    int quenchElectrons(const int itmax, double& etot);
    void dumpRestart();
    void updateRefs();
    void setQuenchTol() const;
    void updatePotAndMasks();
    void setForces(const std::vector<std::vector<double>>& f);
    bool lbfgsLastStepNotAccepted() const;
    void addImageInteractionEnergy(const double eii);
};

#endif
