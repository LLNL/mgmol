// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_MVP_DMStrategy_H
#define MGMOL_MVP_DMStrategy_H

#include "DMStrategy.h"
#include "Energy.h"
#include "MGmol.h"
#include "Rho.h"

#include <iostream>
#include <mpi.h>
#include <vector>

class ProjectedMatricesInterface;
class Ions;
class Electrostatic;

template <class T>
class MVP_DMStrategy : public DMStrategy
{
private:
    T* orbitals_;
    ProjectedMatricesInterface* proj_matrices_;

    MPI_Comm comm_;
    std::ostream& os_;

    Ions& ions_;
    Rho<T>* rho_;
    Energy<T>* energy_;
    Electrostatic* electrostat_;
    const std::vector<std::vector<int>>& global_indexes_;
    MGmol<T>* mgmol_strategy_;

    bool use_old_dm_;

public:
    MVP_DMStrategy(MPI_Comm comm, std::ostream& os, Ions& ions, Rho<T>* rho,
        Energy<T>* energy, Electrostatic* electrostat, MGmol<T>* mgmol_strategy,
        T* orbitals, ProjectedMatricesInterface* proj_matrices,
        const bool use_old_dm);

    void initialize(){};
    int update();

    // H is updated with MVP loop, so no need to compute it outside
    bool needH() const { return false; }

    void stripDM();

    void dressDM();

    void reset() {}
};

#endif
