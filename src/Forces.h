// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_FORCES_H
#define MGMOL_FORCES_H

#include "Hamiltonian.h"
#include "Rho.h"
#include "global.h"

#include <array>
#include <functional>
#include <vector>

#define MGMOL_FD_NPTS 2
#define MGMOL_FD_DELTA 0.002

class Ion;
class Ions;
class ProjectedMatricesInterface;

template <class OrbitalsType>
class Forces
{
private:
    double shift_R_[3 * MGMOL_FD_NPTS][3];

    Hamiltonian<OrbitalsType>* hamiltonian_;
    Rho<OrbitalsType>* rho_;
    ProjectedMatricesInterface* proj_matrices_;

    static Timer lforce_tm_;
    static Timer nlforce_tm_;
    static Timer evaluateShiftedFields_tm_;
    static Timer computeIntegrals_tm_;
    static Timer consolidate_data_;
    static Timer lforce_local_tm_;
    static Timer kbpsi_tm_;
    static Timer energy_tm_;
    static Timer total_tm_;

    void integrals_ion(const Ion& ion, const std::vector<RHODTYPE>& rho,
        const std::vector<POTDTYPE>& vh_rho,
        std::array<double, 3 * MGMOL_FD_NPTS>& integrals,
        const char flag_filter);
    void computeIntegrals(const std::vector<double>& field,
        const std::vector<std::array<double, 3 * MGMOL_FD_NPTS>>& shifted_field,
        std::array<double, 3 * MGMOL_FD_NPTS>& integrals);
    void evaluateShiftedFields(const Ion& ion,
        std::vector<std::array<double, 3 * MGMOL_FD_NPTS>>& var_pot,
        std::vector<std::array<double, 3 * MGMOL_FD_NPTS>>& rhos,
        const char flag_filter);
    void evaluateSupersampledRadialFunc(const std::vector<Vector3D>& positions,
        const double lrad,
        std::vector<std::array<double, 3 * MGMOL_FD_NPTS>>& var,
        std::function<double(double)> const&);
    void evaluateRadialFunc(const std::vector<Vector3D>& positions,
        const double lrad,
        std::vector<std::array<double, 3 * MGMOL_FD_NPTS>>& var,
        std::function<double(double)> const&);
    SquareLocalMatrices<double, MemorySpace::Host> getReplicatedDM();

    void nlforceSparse(OrbitalsType& orbitals, Ions& ions);
    void lforce(Ions& ions, const std::vector<RHODTYPE>& rho,
        const std::vector<POTDTYPE>& vh_rho);
    void external_force(Ions& ions);
    void efield_force(Ions& ions);

public:
    Forces(Hamiltonian<OrbitalsType>* hamiltonian, Rho<OrbitalsType>* rho,
        ProjectedMatricesInterface* proj_matrices);

    void force(OrbitalsType& orbitals, Ions& ions);

    void printTimers(std::ostream& os)
    {
        lforce_tm_.print(os);
        nlforce_tm_.print(os);
        evaluateShiftedFields_tm_.print(os);
        computeIntegrals_tm_.print(os);
        consolidate_data_.print(os);
        lforce_local_tm_.print(os);
        kbpsi_tm_.print(os);
        energy_tm_.print(os);
        total_tm_.print(os);
    }
};

template <class OrbitalsType>
Timer Forces<OrbitalsType>::lforce_tm_("Forces::lforce");
template <class OrbitalsType>
Timer Forces<OrbitalsType>::nlforce_tm_("Forces::nlforce");
template <class OrbitalsType>
Timer Forces<OrbitalsType>::evaluateShiftedFields_tm_(
    "Forces::evalShiftedFields");
template <class OrbitalsType>
Timer Forces<OrbitalsType>::computeIntegrals_tm_("Forces::computeIntegrals");
template <class OrbitalsType>
Timer Forces<OrbitalsType>::consolidate_data_("Forces::consolidate");
template <class OrbitalsType>
Timer Forces<OrbitalsType>::lforce_local_tm_("Forces::lforce_local");
template <class OrbitalsType>
Timer Forces<OrbitalsType>::total_tm_("Forces::total");
template <class OrbitalsType>
Timer Forces<OrbitalsType>::kbpsi_tm_("Forces::KBpsi");
template <class OrbitalsType>
Timer Forces<OrbitalsType>::energy_tm_("Forces::nl_energy");

#endif
