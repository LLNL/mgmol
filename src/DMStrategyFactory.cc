#include "DMStrategyFactory.h"
#include "ReplicatedMatrix.h"

#ifdef MGMOL_USE_SCALAPACK
template <>
DMStrategy<LocGridOrbitals<ORBDTYPE>>*
DMStrategyFactory<LocGridOrbitals<ORBDTYPE>,
    dist_matrix::DistMatrix<double>>::createHamiltonianMVP_DMStrategy(MPI_Comm
                                                                          comm,
    std::ostream& os, Ions& ions, Rho<LocGridOrbitals<ORBDTYPE>>* rho,
    Energy<LocGridOrbitals<ORBDTYPE>>* energy, Electrostatic* electrostat,
    Hamiltonian<LocGridOrbitals<ORBDTYPE>>* hamiltonian,
    MGmol<LocGridOrbitals<ORBDTYPE>>* mgmol_strategy,
    ProjectedMatricesInterface* /*proj_matrices*/,
    LocGridOrbitals<ORBDTYPE>* orbitals, const bool short_sighted)
{
    if (short_sighted)
    {
        DMStrategy<LocGridOrbitals<ORBDTYPE>>* dm_strategy
            = new HamiltonianMVP_DMStrategy<VariableSizeMatrix<sparserow>,
                ProjectedMatricesSparse, LocGridOrbitals<ORBDTYPE>>(comm, os,
                ions, rho, energy, electrostat, hamiltonian, mgmol_strategy,
                orbitals);

        return dm_strategy;
    }
    else
    {
        DMStrategy<LocGridOrbitals<ORBDTYPE>>* dm_strategy
            = new HamiltonianMVP_DMStrategy<
                dist_matrix::DistMatrix<DISTMATDTYPE>,
                ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>,
                LocGridOrbitals<ORBDTYPE>>(comm, os, ions, rho, energy,
                electrostat, hamiltonian, mgmol_strategy, orbitals);

        return dm_strategy;
    }
}
#endif

template <>
DMStrategy<LocGridOrbitals<ORBDTYPE>>*
DMStrategyFactory<LocGridOrbitals<ORBDTYPE>,
    ReplicatedMatrix>::createHamiltonianMVP_DMStrategy(MPI_Comm comm,
    std::ostream& /*os*/, Ions& /*ions*/,
    Rho<LocGridOrbitals<ORBDTYPE>>* /*rho*/,
    Energy<LocGridOrbitals<ORBDTYPE>>* /*energy*/,
    Electrostatic* /*electrostat*/,
    Hamiltonian<LocGridOrbitals<ORBDTYPE>>* /*hamiltonian*/,
    MGmol<LocGridOrbitals<ORBDTYPE>>* /*mgmol_strategy*/,
    ProjectedMatricesInterface* /*proj_matrices*/,
    LocGridOrbitals<ORBDTYPE>* /*orbitals*/, const bool /*short_sighted*/)
{

    std::cerr << "DMStrategy not implemented" << std::endl;
    MPI_Abort(comm, EXIT_FAILURE);

    return nullptr;
}

#ifdef MGMOL_USE_SCALAPACK
template <>
DMStrategy<ExtendedGridOrbitals<ORBDTYPE>>*
DMStrategyFactory<ExtendedGridOrbitals<ORBDTYPE>,
    dist_matrix::DistMatrix<double>>::createHamiltonianMVP_DMStrategy(MPI_Comm
                                                                          comm,
    std::ostream& os, Ions& ions, Rho<ExtendedGridOrbitals<ORBDTYPE>>* rho,
    Energy<ExtendedGridOrbitals<ORBDTYPE>>* energy, Electrostatic* electrostat,
    Hamiltonian<ExtendedGridOrbitals<ORBDTYPE>>* hamiltonian,
    MGmol<ExtendedGridOrbitals<ORBDTYPE>>* mgmol_strategy,
    ProjectedMatricesInterface* /*proj_matrices*/,
    ExtendedGridOrbitals<ORBDTYPE>* orbitals, const bool short_sighted)
{
    (void)short_sighted;

    DMStrategy<ExtendedGridOrbitals<ORBDTYPE>>* dm_strategy
        = new HamiltonianMVP_DMStrategy<dist_matrix::DistMatrix<DISTMATDTYPE>,
            ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>,
            ExtendedGridOrbitals<ORBDTYPE>>(comm, os, ions, rho, energy,
            electrostat, hamiltonian, mgmol_strategy, orbitals);

    return dm_strategy;
}
#endif

template <>
DMStrategy<ExtendedGridOrbitals<ORBDTYPE>>*
DMStrategyFactory<ExtendedGridOrbitals<ORBDTYPE>,
    ReplicatedMatrix>::createHamiltonianMVP_DMStrategy(MPI_Comm comm,
    std::ostream& os, Ions& ions, Rho<ExtendedGridOrbitals<ORBDTYPE>>* rho,
    Energy<ExtendedGridOrbitals<ORBDTYPE>>* energy, Electrostatic* electrostat,
    Hamiltonian<ExtendedGridOrbitals<ORBDTYPE>>* hamiltonian,
    MGmol<ExtendedGridOrbitals<ORBDTYPE>>* mgmol_strategy,
    ProjectedMatricesInterface* /*proj_matrices*/,
    ExtendedGridOrbitals<ORBDTYPE>* orbitals, const bool short_sighted)
{
    (void)short_sighted;

    DMStrategy<ExtendedGridOrbitals<ORBDTYPE>>* dm_strategy
        = new HamiltonianMVP_DMStrategy<ReplicatedMatrix,
            ProjectedMatrices<ReplicatedMatrix>,
            ExtendedGridOrbitals<ORBDTYPE>>(comm, os, ions, rho, energy,
            electrostat, hamiltonian, mgmol_strategy, orbitals);

    return dm_strategy;
}
