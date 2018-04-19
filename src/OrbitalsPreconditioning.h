// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef OrbitalsPreconditioning_H
#define OrbitalsPreconditioning_H

#include "GridFuncVector.h"
#include "Lap.h"
#include "LocGridOrbitals.h"
#include "Preconditioning.h"

class Masks4Orbitals;
//class LocGridOrbitals;
class MasksSet;
class ProjectedMatricesInterface;
class Potentials;
class LocalizationRegions;

class OrbitalsPreconditioning
{
private:
    Preconditioning<MGPRECONDTYPE>* precond_;
    pb::GridFuncVector<MGPRECONDTYPE>* gfv_work_;
    
    pb::GridFuncVector<MGPRECONDTYPE>* data_wghosts_;

    // coefficient for preconditioning
    double  gamma_;
    
    bool is_set_;
    
    //preconditioner precision != orbitals precision
    bool mixed_precision_;

    // timers
    static Timer   precond_tm_;    

public:
    OrbitalsPreconditioning()
    {
        is_set_=false;
        precond_=0;
        gfv_work_=0;
        data_wghosts_=0;
        
        mixed_precision_= (sizeof(MGPRECONDTYPE)!=sizeof(ORBDTYPE));
    };
    
    ~OrbitalsPreconditioning();
    
    void setup(LocGridOrbitals& orbitals,
               const short mg_levels,
               const short lap_type,
               MasksSet*,
               LocalizationRegions*);
    void precond_mg(LocGridOrbitals& orbitals);
    void setGamma(const pb::Lap<ORBDTYPE>& lapOper, const Potentials& pot,
                  const short mg_levels,
                  ProjectedMatricesInterface* proj_matrices);
    static void printTimers(std::ostream& os);
};

#endif
