// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_IonicAlgorithm_H
#define MGMOL_IonicAlgorithm_H

#include "ConstraintSet.h"
#include "IonicStepper.h"
#include "Ions.h"
#include "LocalizationRegions.h"
#include "MGmol.h"
#include "Rho.h"

#include <vector>

class MasksSet;

template <class OrbitalsType>
class IonicAlgorithm
{
private:
    OrbitalsType** orbitals_;
    Ions& ions_;
    Rho<OrbitalsType>& rho_;
    ConstraintSet& constraints_;
    std::shared_ptr<LocalizationRegions> lrs_;
    MasksSet& masks_;

    MGmol<OrbitalsType>& mgmol_strategy_;
    IonicStepper* stepper_;
    std::vector<std::string>& ions_names_;

protected:
    std::vector<double>& tau0_;
    std::vector<double>& taup_;
    std::vector<double>& fion_;
    std::vector<double>& pmass_;
    const std::vector<short>& atmove_;
    const std::vector<int>& gid_;

    void registerStepper(IonicStepper* stepper);

public:
    IonicAlgorithm(OrbitalsType** orbitals, Ions& ions, Rho<OrbitalsType>& rho,
        ConstraintSet& constraints, std::shared_ptr<LocalizationRegions> lrs,
        MasksSet& masks, MGmol<OrbitalsType>&);

    virtual ~IonicAlgorithm() {};

    virtual void init(HDFrestart* h5f_file);
    virtual int run1step();
    void setupConstraints();
    void computeForces();
    bool checkTolForces(const double tol);
    void dumpRestart();
    void updatePotAndMasks();
    void setForces(const std::vector<std::vector<double>>& f);
    virtual int quenchElectrons(const int itmax, double& etot);
};

#endif
