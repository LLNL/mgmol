// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_FORCES_H
#define MGMOL_FORCES_H

#include "Hamiltonian.h"
#include "Rho.h"
#include "global.h"

#define NPTS 2
#define DELTAC 0.002

class Ion;
class Ions;
class LocGridOrbitals;
class ProjectedMatricesInterface;

class Forces
{
private:
    Hamiltonian<LocGridOrbitals>* hamiltonian_;
    Rho<LocGridOrbitals>* rho_;
    ProjectedMatricesInterface* proj_matrices_;

    static Timer lforce_tm_;
    static Timer nlforce_tm_;
    static Timer get_var_tm_;
    static Timer get_loc_proj_tm_;
    static Timer consolidate_data_;
    static Timer lforce_local_tm_;
    static Timer kbpsi_tm_;
    static Timer energy_tm_;
    static Timer total_tm_;

    void lforce_ion(Ion& ion, RHODTYPE* rho, double** loc_proj);
    void get_loc_proj(RHODTYPE* rho, const int* const pvec, double*** var_pot,
        double*** var_charge, const int docount, double** loc_proj);
    int get_var(Ion& ion, int* pvec, double*** var_pot, double*** var_charge);

public:
    Forces(Hamiltonian<LocGridOrbitals>* hamiltonian, Rho<LocGridOrbitals>* rho,
        ProjectedMatricesInterface* proj_matrices)
        : hamiltonian_(hamiltonian), rho_(rho), proj_matrices_(proj_matrices)
    {
        assert(hamiltonian_ != 0);
        assert(rho_ != 0);
        assert(proj_matrices_ != 0);
    }

    void nlforce(LocGridOrbitals& orbitals, Ions& ions);
    void nlforceSparse(LocGridOrbitals& orbitals, Ions& ions);
    void lforce(Ions& ions, RHODTYPE* rho);
    void force(LocGridOrbitals& orbitals, Ions& ions);

    void printTimers(ostream& os)
    {
        lforce_tm_.print(os);
        nlforce_tm_.print(os);
        get_var_tm_.print(os);
        get_loc_proj_tm_.print(os);
        consolidate_data_.print(os);
        lforce_local_tm_.print(os);
        kbpsi_tm_.print(os);
        energy_tm_.print(os);
        total_tm_.print(os);
    }
};

#endif
