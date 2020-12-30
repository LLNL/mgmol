// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Hamiltonian.h"
#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"
#include "Potentials.h"
#include "ProjectedMatrices.h"
#include "ReplicatedMatrix.h"

template <class T>
Hamiltonian<T>::Hamiltonian()
{
    itindex_ = -1;
    lapOper_ = nullptr;
    hlphi_   = nullptr;
    pot_     = new Potentials();
}

template <class T>
Hamiltonian<T>::~Hamiltonian()
{
    if (hlphi_ != nullptr) delete hlphi_;
    if (lapOper_ != nullptr) delete lapOper_;
    delete pot_;
}

template <class T>
void Hamiltonian<T>::setup(const pb::Grid& myGrid, const int lap_type)
{
    if (lapOper_ != nullptr) delete lapOper_;
    lapOper_ = LapFactory<ORBDTYPE>::createLap(myGrid, lap_type);
}

template <class T>
const T& Hamiltonian<T>::applyLocal(T& phi, const bool force)
{
    assert(phi.getIterativeIndex() >= 0);
    assert(pot_->getIterativeIndex() >= 0);

    if (hlphi_ == nullptr) hlphi_ = new T("Hphi", phi, false);
    if (!hlphi_->isCompatibleWith(phi))
    {
        delete hlphi_;
        itindex_ = -1;
        hlphi_   = new T("Hphi", phi, false);
    }
    const int new_index
        = 100 * phi.getIterativeIndex() + pot_->getIterativeIndex();
#ifdef DEBUG
    if (onpe0)
    {
        (*MPIdata::sout) << "Hamiltonian<T>::applyLocal(), new_index ="
                         << new_index << endl;
        (*MPIdata::sout) << "Hamiltonian<T>::applyLocal(), itindex_  ="
                         << itindex_ << endl;
    }
#endif
    if (force || new_index != itindex_)
    {
        applyLocal(phi.chromatic_number(), phi, *hlphi_);

        itindex_ = new_index;
#ifdef PRINT_OPERATIONS
    }
    else
    {
        if (onpe0)
            (*MPIdata::sout)
                << "Hamiltonian<T>::hlphi up to date, itindex_=" << itindex_
                << endl;
#endif
    }
    return *hlphi_;
}

template <class T>
void Hamiltonian<T>::applyLocal(const int ncolors, T& phi, T& hphi)
{
    apply_Hloc_tm_.start();
#ifdef PRINT_OPERATIONS
    if (onpe0)
        (*MPIdata::sout) << "Hamiltonian<T>::applyLocal() for " << ncolors
                         << " states" << endl;
#endif

    const Control& ct      = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    const POTDTYPE* const vtot = pot_->vtot();

    phi.setDataWithGhosts();
    phi.trade_boundaries();

    using memory_space_type = typename T::memory_space_type;

    if (ct.Mehrstellen())
    {
        pb::GridFunc<POTDTYPE> gfpot(
            mygrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
        gfpot.assign(vtot);
        if (ct.Mehrstellen()) gfpot.trade_boundaries();
        const std::vector<std::vector<int>>& gid(phi.getOverlappingGids());
        pb::GridFuncVector<ORBDTYPE, memory_space_type> gfvw1(
            mygrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2], gid);
        pb::GridFuncVector<ORBDTYPE, memory_space_type>& gfvphi(
            *phi.getPtDataWGhosts());
        gfvw1.prod(gfvphi, gfpot);

        pb::GridFunc<ORBDTYPE> gf_work1(
            mygrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
        pb::GridFunc<ORBDTYPE> gf_work2(
            mygrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
        for (int i = 0; i < ncolors; i++)
        {
            // work1 = B*V*psi
            lapOper_->rhs(gfvw1.getGridFunc(i), gf_work1);

            // work2 = -Lap*phi
            lapOper_->apply(phi.getFuncWithGhosts(i), gf_work2);

            gf_work1 += gf_work2;
            hphi.setPsi(gf_work1, i);
        }
    }
    else
    {
#pragma omp parallel for
        for (int i = 0; i < ncolors; i++)
        {
            using memory_space_type   = typename T::memory_space_type;
            ORBDTYPE* ihphi           = hphi.getPsi(i);
            unsigned int const size   = hphi.getNumpt();
            ORBDTYPE* ihphi_host_view = MemorySpace::Memory<ORBDTYPE,
                memory_space_type>::allocate_host_view(size);
            MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
                hphi.getPsi(i), size, ihphi_host_view);

            lapOper_->applyWithPot(
                phi.getFuncWithGhosts(i), vtot, ihphi_host_view);

            MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_dev(
                ihphi_host_view, size, ihphi);
            MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
                ihphi_host_view);
        }
    }

    apply_Hloc_tm_.stop();
}

// add to hij the elements <phi1|Hloc|phi2>
// corresponding to the local part of the Hamiltonian
template <>
template <>
void Hamiltonian<LocGridOrbitals>::addHlocal2matrix(LocGridOrbitals& phi1,
    LocGridOrbitals& phi2, dist_matrix::DistMatrix<double>& hij,
    const bool force)
{
    applyLocal(phi2, force);

#ifdef PRINT_OPERATIONS
    if (onpe0) (*MPIdata::sout) << "Hamiltonian<T>::addHlocal2matrix()" << endl;
#endif

    phi1.addDotWithNcol2Matrix(*hlphi_, hij);
}

template <>
template <>
void Hamiltonian<ExtendedGridOrbitals>::addHlocal2matrix(
    ExtendedGridOrbitals& phi1, ExtendedGridOrbitals& phi2,
    dist_matrix::DistMatrix<double>& hij, const bool force)
{
    applyLocal(phi2, force);

#ifdef PRINT_OPERATIONS
    if (onpe0) (*MPIdata::sout) << "Hamiltonian<T>::addHlocal2matrix()" << endl;
#endif

    // hij.print(std::cout, 0, 0, 5, 5);

    phi1.addDotWithNcol2Matrix(*hlphi_, hij);

    // hij.print(std::cout, 0, 0, 5, 5);
}

#ifdef HAVE_MAGMA
template <>
template <>
void Hamiltonian<ExtendedGridOrbitals>::addHlocal2matrix(
    ExtendedGridOrbitals& phi1, ExtendedGridOrbitals& phi2,
    ReplicatedMatrix& hij, const bool force)
{
    applyLocal(phi2, force);

#ifdef PRINT_OPERATIONS
    if (onpe0) (*MPIdata::sout) << "Hamiltonian<T>::addHlocal2matrix()" << endl;
#endif

    phi1.addDotWithNcol2Matrix(*hlphi_, hij);
}
#endif

template <class T>
void Hamiltonian<T>::addHlocalij(
    T& phi1, T& phi2, ProjectedMatricesInterface* proj_matrices)
{
    applyLocal(phi2);

#ifdef PRINT_OPERATIONS
    if (onpe0) (*MPIdata::sout) << "Hamiltonian<T>::addHLocalij()" << endl;
#endif

    addHlocalij(phi1, proj_matrices);
}

template <class T>
void Hamiltonian<T>::addHlocalij(
    T& phi1, ProjectedMatricesInterface* proj_matrices)
{
    SquareLocalMatrices<MATDTYPE, MemorySpace::Host> slm(
        phi1.subdivx(), phi1.chromatic_number());

    phi1.computeLocalProduct(*hlphi_, slm);

    proj_matrices->setLocalMatrixElementsHl(slm);

    proj_matrices->consolidateH();
}

template <>
template <>
void Hamiltonian<LocGridOrbitals>::addHlocal2matrix(LocGridOrbitals& phi1,
    LocGridOrbitals& phi2, VariableSizeMatrix<sparserow>& mat, const bool force)
{
    Control& ct = *(Control::instance());

    applyLocal(phi2, force);

#ifdef PRINT_OPERATIONS
    if (onpe0) (*MPIdata::sout) << "Hamiltonian<T>::addHLocalij()" << endl;
#endif

    SquareLocalMatrices<MATDTYPE, MemorySpace::Host> ss(
        phi1.subdivx(), phi1.chromatic_number());

    phi1.computeLocalProduct(*hlphi_, ss);

    mat.insertMatrixElements(ss, phi1.getOverlappingGids(), ct.numst);
}

template Hamiltonian<LocGridOrbitals>::Hamiltonian();
template Hamiltonian<ExtendedGridOrbitals>::Hamiltonian();

template Hamiltonian<LocGridOrbitals>::~Hamiltonian();
template Hamiltonian<ExtendedGridOrbitals>::~Hamiltonian();

template void Hamiltonian<LocGridOrbitals>::setup(pb::Grid const&, int);
template void Hamiltonian<ExtendedGridOrbitals>::setup(pb::Grid const&, int);

template const LocGridOrbitals& Hamiltonian<LocGridOrbitals>::applyLocal(
    LocGridOrbitals&, const bool);
template const ExtendedGridOrbitals&
Hamiltonian<ExtendedGridOrbitals>::applyLocal(
    ExtendedGridOrbitals&, const bool);
template void Hamiltonian<LocGridOrbitals>::addHlocalij(LocGridOrbitals&,
    LocGridOrbitals&, ProjectedMatricesInterface* proj_matrices);
template void Hamiltonian<ExtendedGridOrbitals>::addHlocalij(
    ExtendedGridOrbitals&, ExtendedGridOrbitals&,
    ProjectedMatricesInterface* proj_matrices);
template void Hamiltonian<LocGridOrbitals>::addHlocalij(
    LocGridOrbitals&, ProjectedMatricesInterface* proj_matrices);
template void Hamiltonian<ExtendedGridOrbitals>::addHlocalij(
    ExtendedGridOrbitals&, ProjectedMatricesInterface* proj_matrices);
template void Hamiltonian<LocGridOrbitals>::addHlocal2matrix(LocGridOrbitals&,
    LocGridOrbitals&, VariableSizeMatrix<sparserow>& mat, const bool force);
