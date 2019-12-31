#include "DMStrategyFactory.h"

template <>
DMStrategy* DMStrategyFactory<LocGridOrbitals>::createHamiltonianMVP_DMStrategy(
    MPI_Comm comm, std::ostream& os, Ions& ions, Rho<LocGridOrbitals>* rho,
    Energy<LocGridOrbitals>* energy, Electrostatic* electrostat,
    MGmol<LocGridOrbitals>* mgmol_strategy,
    ProjectedMatricesInterface* proj_matrices, LocGridOrbitals* orbitals,
    const bool short_sighted)
{
    if (short_sighted)
    {
        DMStrategy* dm_strategy
            = new HamiltonianMVP_DMStrategy<VariableSizeMatrix<sparserow>,
                ProjectedMatricesSparse, LocGridOrbitals>(comm, os, ions, rho,
                energy, electrostat, mgmol_strategy, orbitals);

        return dm_strategy;
    }
    else
    {
        DMStrategy* dm_strategy = new HamiltonianMVP_DMStrategy<
            dist_matrix::DistMatrix<DISTMATDTYPE>, ProjectedMatrices,
            LocGridOrbitals>(
            comm, os, ions, rho, energy, electrostat, mgmol_strategy, orbitals);

        return dm_strategy;
    }
}

template <>
DMStrategy*
DMStrategyFactory<ExtendedGridOrbitals>::createHamiltonianMVP_DMStrategy(
    MPI_Comm comm, std::ostream& os, Ions& ions, Rho<ExtendedGridOrbitals>* rho,
    Energy<ExtendedGridOrbitals>* energy, Electrostatic* electrostat,
    MGmol<ExtendedGridOrbitals>* mgmol_strategy,
    ProjectedMatricesInterface* proj_matrices, ExtendedGridOrbitals* orbitals,
    const bool short_sighted)
{
    (void)short_sighted;

    DMStrategy* dm_strategy
        = new HamiltonianMVP_DMStrategy<dist_matrix::DistMatrix<DISTMATDTYPE>,
            ProjectedMatrices, ExtendedGridOrbitals>(
            comm, os, ions, rho, energy, electrostat, mgmol_strategy, orbitals);

    return dm_strategy;
}
