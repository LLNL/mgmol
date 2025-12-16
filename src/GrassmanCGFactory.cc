#include "GrassmanCGFactory.h"
#include "GrassmanCG.h"
#include "GrassmanCGSparse.h"
#include "LocGridOrbitals.h"

template <>
OrbitalsStepper<LocGridOrbitals<ORBDTYPE>>*
GrassmanCGFactory<LocGridOrbitals<ORBDTYPE>>::create(
    Hamiltonian<LocGridOrbitals<ORBDTYPE>>* hamiltonian,
    ProjectedMatricesInterface* proj_matrices,
    MGmol<LocGridOrbitals<ORBDTYPE>>* mgmol_strategy, Ions& ions,
    std::ostream& os, const bool short_sighted)
{
    OrbitalsStepper<LocGridOrbitals<ORBDTYPE>>* stepper;

    if (short_sighted)
    {
        stepper = new GrassmanCGSparse<LocGridOrbitals<ORBDTYPE>>(
            hamiltonian, proj_matrices, mgmol_strategy, ions, os);
    }
    else
    {
        stepper = new GrassmanCG<LocGridOrbitals<ORBDTYPE>>(
            hamiltonian, proj_matrices, mgmol_strategy, ions, os);
    }

    return stepper;
}

template <>
OrbitalsStepper<ExtendedGridOrbitals<ORBDTYPE>>*
GrassmanCGFactory<ExtendedGridOrbitals<ORBDTYPE>>::create(
    Hamiltonian<ExtendedGridOrbitals<ORBDTYPE>>* hamiltonian,
    ProjectedMatricesInterface* proj_matrices,
    MGmol<ExtendedGridOrbitals<ORBDTYPE>>* mgmol_strategy, Ions& ions,
    std::ostream& os, const bool /*short_sighted*/)
{
    OrbitalsStepper<ExtendedGridOrbitals<ORBDTYPE>>* stepper
        = new GrassmanCG<ExtendedGridOrbitals<ORBDTYPE>>(
            hamiltonian, proj_matrices, mgmol_strategy, ions, os);

    return stepper;
}
