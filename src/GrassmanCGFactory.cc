#include "GrassmanCGFactory.h"
#include "GrassmanCGSparse.h"
#include "GrassmanCG.h"
#include "LocGridOrbitals.h"

template <>
OrbitalsStepper<LocGridOrbitals>* GrassmanCGFactory<LocGridOrbitals>::create(
    Hamiltonian<LocGridOrbitals>* hamiltonian,
    ProjectedMatricesInterface* proj_matrices,
    MGmol<LocGridOrbitals>* mgmol_strategy,
    Ions& ions, std::ostream& os, const bool short_sighted)
{
    OrbitalsStepper<LocGridOrbitals>* stepper;

    if (short_sighted)
    {
        stepper = new GrassmanCGSparse<LocGridOrbitals>(
            hamiltonian, proj_matrices, mgmol_strategy, ions, os);
    }
    else
    {
        stepper = new GrassmanCG<LocGridOrbitals>(
            hamiltonian, proj_matrices, mgmol_strategy, ions, os);
    }

    return stepper;
}

template <>
OrbitalsStepper<ExtendedGridOrbitals>*
GrassmanCGFactory<ExtendedGridOrbitals>::create(
    Hamiltonian<ExtendedGridOrbitals>* hamiltonian,
    ProjectedMatricesInterface* proj_matrices,
    MGmol<ExtendedGridOrbitals>* mgmol_strategy,
    Ions& ions, std::ostream& os, const bool short_sighted)
{
    OrbitalsStepper<ExtendedGridOrbitals>* stepper
        = new GrassmanCG<ExtendedGridOrbitals>(
            hamiltonian, proj_matrices, mgmol_strategy, ions, os);

    return stepper;
}

