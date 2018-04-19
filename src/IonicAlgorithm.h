// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#ifndef IonicAlgorithm_H
#define IonicAlgorithm_H

#include "IonicStepper.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "Rho.h"
#include "ConstraintSet.h"
#include "LocalizationRegions.h"

#include <vector>

class MasksSet;
class Electrostatic;
class MGmol;

class IonicAlgorithm
{
private:
    LocGridOrbitals** orbitals_;
    Ions& ions_;
    Rho& rho_;
    ConstraintSet& constraints_;
    LocalizationRegions& lrs_;
    MasksSet& masks_;
    
    MGmol& mgmol_strategy_;
    IonicStepper* stepper_;
    std::vector<std::string>& ions_names_;

protected:
    std::vector<double>& tau0_;   
    std::vector<double>& taup_;   
    std::vector<double>& fion_;   
    std::vector<double>& pmass_;  
    const std::vector<short>&  atmove_; 
    const std::vector<int>& gid_;
    
    void registerStepper(IonicStepper* stepper);
    
public:
    IonicAlgorithm(LocGridOrbitals** orbitals,
          Ions& ions,
          Rho& rho,
          ConstraintSet& constraints,
          LocalizationRegions& lrs,
          MasksSet& masks,
          MGmol&);
    
    virtual ~IonicAlgorithm(){};
    
    virtual void init(HDFrestart* h5f_file);
    virtual int run1step();
    void setupConstraints();
    void computeForces();
    bool checkTolForces(const double tol);
    void dumpRestart();
    void updatePotAndMasks();
    void setForces(const std::vector<std::vector<double> >& f);
    virtual int quenchElectrons(const int itmax, double& etot);
};

#endif
