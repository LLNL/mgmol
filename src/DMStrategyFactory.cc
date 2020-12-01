#include "DMStrategyFactory.h"
#include "ReplicatedMatrix.h"

template <>
DMStrategy* DMStrategyFactory<LocGridOrbitals,dist_matrix::DistMatrix<double>>::
    createHamiltonianMVP_DMStrategy(
    MPI_Comm comm, std::ostream& os, Ions& ions, Rho<LocGridOrbitals>* rho,
    Energy<LocGridOrbitals>* energy, Electrostatic* electrostat,
    MGmol<LocGridOrbitals>* mgmol_strategy,
    ProjectedMatricesInterface* /*proj_matrices*/, LocGridOrbitals* orbitals,
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
            dist_matrix::DistMatrix<DISTMATDTYPE>,
            ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>,
            LocGridOrbitals>(
            comm, os, ions, rho, energy, electrostat, mgmol_strategy, orbitals);

        return dm_strategy;
    }
}

template <>
DMStrategy*
DMStrategyFactory<ExtendedGridOrbitals, dist_matrix::DistMatrix<double>>::createHamiltonianMVP_DMStrategy(
    MPI_Comm comm, std::ostream& os, Ions& ions, Rho<ExtendedGridOrbitals>* rho,
    Energy<ExtendedGridOrbitals>* energy, Electrostatic* electrostat,
    MGmol<ExtendedGridOrbitals>* mgmol_strategy,
    ProjectedMatricesInterface* /*proj_matrices*/,
    ExtendedGridOrbitals* orbitals, const bool short_sighted)
{
    (void)short_sighted;

    DMStrategy* dm_strategy
        = new HamiltonianMVP_DMStrategy<dist_matrix::DistMatrix<DISTMATDTYPE>,
            ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>,
            ExtendedGridOrbitals>(
            comm, os, ions, rho, energy, electrostat, mgmol_strategy, orbitals);

    return dm_strategy;
}

#ifdef HAVE_MAGMA
template <>
DMStrategy*
DMStrategyFactory<ExtendedGridOrbitals, ReplicatedMatrix>::createHamiltonianMVP_DMStrategy(
    MPI_Comm comm, std::ostream& os, Ions& ions, Rho<ExtendedGridOrbitals>* rho,
    Energy<ExtendedGridOrbitals>* energy, Electrostatic* electrostat,
    MGmol<ExtendedGridOrbitals>* mgmol_strategy,
    ProjectedMatricesInterface* /*proj_matrices*/,
    ExtendedGridOrbitals* orbitals, const bool short_sighted)
{
    (void)short_sighted;

    DMStrategy* dm_strategy
        = new HamiltonianMVP_DMStrategy<ReplicatedMatrix,
            ProjectedMatrices<ReplicatedMatrix>,
            ExtendedGridOrbitals>(
            comm, os, ions, rho, energy, electrostat, mgmol_strategy, orbitals);

    return dm_strategy;
}
#endif
