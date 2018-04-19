// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef HamiltonianMVP_DMStrategy_H
#define HamiltonianMVP_DMStrategy_H

#include "DMStrategy.h"
#include "HamiltonianMVPSolver.h"
class ProjectedMatricesInterface;
class LocGridOrbitals;
class Ions;
class Hamiltonian;
class Rho;
class Energy;
class MGmol;
class Electrostatic;
class MGmol;
//class HamiltonianMVPSolver;

template <class T1, class T2, class T3>
class HamiltonianMVP_DMStrategy:public DMStrategy
{
private:
    LocGridOrbitals* orbitals_;

    MPI_Comm comm_;
    std::ostream& os_;

    Ions& ions_;
    Rho* rho_;
    Energy* energy_;
    Electrostatic* electrostat_;
    const std::vector<std::vector<int> >& global_indexes_;
    MGmol* mgmol_strategy_;
    
    HamiltonianMVPSolver<T1, T2, T3>* solver_;
    
public:
    HamiltonianMVP_DMStrategy(MPI_Comm comm, std::ostream& os, 
        Ions& ions,
        Rho* rho,
        Energy* energy,
        Electrostatic* electrostat,
        MGmol* mgmol_strategy,
        LocGridOrbitals* orbitals);
    
    ~HamiltonianMVP_DMStrategy();

    void initialize();
    int update();

    //H is updated with HamiltonianMVP loop, so no need to compute it outside
    bool needH()const{ return false; }

    void stripDM();
    
    void dressDM();
    
    void reset();
};

#endif
