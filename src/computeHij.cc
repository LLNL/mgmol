// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// OrbitalsType& This file is part of MGmol. For details, see
// https://github.com/llnl/mgmol. Please also read this link
// https://github.com/llnl/mgmol/LICENSE

#include <cassert>

#include "Control.h"
#include "Energy.h"
#include "GridFuncVector.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "KBPsiMatrixSparse.h"
#include "Lap.h"
#include "MGmol.h"
#include "MGmol_MPI.h"
#include "MGmol_blas1.h"
#include "Mesh.h"
#include "Potentials.h"
#include "ProjectedMatricesInterface.h"
#include "ProjectedMatricesSparse.h"
#include "ReplicatedMatrix.h"

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrix.h"
#include "SquareSubMatrix2DistMatrix.h"
#endif

template <>
template <>
void MGmol<LocGridOrbitals<ORBDTYPE>>::computeHij(
    LocGridOrbitals<ORBDTYPE>& orbitals_i,
    LocGridOrbitals<ORBDTYPE>& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi_i,
    const KBPsiMatrixSparse* const kbpsi_j, VariableSizeMatrix<sparserow>& mat,
    const bool consolidate)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHij() at line " << __LINE__ << std::endl;
#endif

    // compute phi_i^T*Hnl*Phi_j
    kbpsi_i->computeHvnlMatrix(kbpsi_j, ions, mat);

    // add local Hamiltonian part to phi_i^T*H*phi_j
    hamiltonian_->addHlocal2matrix(orbitals_i, orbitals_j, mat, true);

    // sum matrix elements among processors
    if (consolidate)
    {
        // get indexes of centered data
        std::vector<int> locfcns;
        (*lrs_).getLocalSubdomainIndices(locfcns);

        // gather/ distribute data from neighbors whose centered functions
        // overlap with functions centered on local subdomain
        Mesh* mymesh             = Mesh::instance();
        const pb::Grid& mygrid   = mymesh->grid();
        const pb::PEenv& myPEenv = mymesh->peenv();
        double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };

        DataDistribution distributorH(
            "prodH", 3. * (*lrs_).max_radii(), myPEenv, domain);
        distributorH.augmentLocalData(mat, false);

        DataDistribution distributor(
            "Hij", 2 * (*lrs_).max_radii(), myPEenv, domain);

        distributor.consolidateMatrix(locfcns, mat);
    }
}

template <>
void MGmol<LocGridOrbitals<ORBDTYPE>>::computeHij(
    LocGridOrbitals<ORBDTYPE>& orbitals_i,
    LocGridOrbitals<ORBDTYPE>& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi, VariableSizeMatrix<sparserow>& mat,
    const bool consolidate)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHij() at line " << __LINE__ << std::endl;
#endif

    kbpsi->computeHvnlMatrix(ions, mat);

    // add local Hamiltonian part to phi^T*H*phi
    hamiltonian_->addHlocal2matrix(orbitals_i, orbitals_j, mat, true);

    // sum matrix elements among processors
    if (consolidate)
    {
        // gather/ distribute data from neighbors whose centered functions
        // overlap with functions centered on local subdomain
        Mesh* mymesh             = Mesh::instance();
        const pb::Grid& mygrid   = mymesh->grid();
        const pb::PEenv& myPEenv = mymesh->peenv();
        double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };

        // get data distribution object
        DataDistribution distributorH(
            "prodH", 3. * (*lrs_).max_radii(), myPEenv, domain);
        // gather data for only centered data
        distributorH.augmentLocalData(mat, false);

        // get indexes of centered data
        std::vector<int> locfcns;
        (*lrs_).getLocalSubdomainIndices(locfcns);

        // sparsify matrix rows according to pattern/ clear rows corresponding
        // to non-centered data
        std::vector<bool> pattern(mat.n(), 0);
        for (std::vector<int>::iterator it = locfcns.begin();
             it != locfcns.end(); ++it)
        {
            const int* rindex = (int*)mat.getTableValue(*it);
            pattern[*rindex]  = 1;
        }
        mat.sparsify(pattern);

        DataDistribution distributor(
            "Hij", 2 * (*lrs_).max_radii(), myPEenv, domain);
        distributor.updateLocalRows(mat, true);
    }
}

#ifdef MGMOL_USE_SCALAPACK
template <class OrbitalsType>
void MGmol<OrbitalsType>::computeHij_private(OrbitalsType& orbitals_i,
    OrbitalsType& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi_i,
    const KBPsiMatrixSparse* const kbpsi_j,
    dist_matrix::DistMatrix<DISTMATDTYPE>& hij)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHij() at line " << __LINE__ << std::endl;
#endif

    hij.clear();

    SquareSubMatrix<double> submat(kbpsi_i->computeHvnlMatrix(kbpsi_j, ions));

    SquareSubMatrix2DistMatrix* ss2dm = SquareSubMatrix2DistMatrix::instance();
    ss2dm->accumulate(submat, hij, 0.);

    // add local Hamiltonian part to phi^T*H*phi
    hamiltonian_->addHlocal2matrix(orbitals_i, orbitals_j, hij, true);
}

template <>
template <>
void MGmol<LocGridOrbitals<ORBDTYPE>>::computeHij(
    LocGridOrbitals<ORBDTYPE>& orbitals_i,
    LocGridOrbitals<ORBDTYPE>& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi,
    const KBPsiMatrixSparse* const kbpsi_j,
    dist_matrix::DistMatrix<DISTMATDTYPE>& hij, const bool consolidate)
{
    (void)consolidate;

    computeHij_private(orbitals_i, orbitals_j, ions, kbpsi, kbpsi_j, hij);
}

template <>
template <>
void MGmol<ExtendedGridOrbitals<ORBDTYPE>>::computeHij(
    ExtendedGridOrbitals<ORBDTYPE>& orbitals_i,
    ExtendedGridOrbitals<ORBDTYPE>& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi,
    const KBPsiMatrixSparse* const kbpsi_j,
    dist_matrix::DistMatrix<DISTMATDTYPE>& hij, const bool consolidate)
{
    (void)consolidate;

    computeHij_private(orbitals_i, orbitals_j, ions, kbpsi, kbpsi_j, hij);
}

template <class OrbitalsType>
void MGmol<OrbitalsType>::computeHij(OrbitalsType& orbitals_i,
    OrbitalsType& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi, dist_matrix::DistMatrix<double>& hij,
    const bool consolidate)
{
    (void)consolidate;

    computeHij_private(orbitals_i, orbitals_j, ions, kbpsi, hij);
}

template <class OrbitalsType>
void MGmol<OrbitalsType>::computeHij_private(OrbitalsType& orbitals_i,
    OrbitalsType& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi, dist_matrix::DistMatrix<double>& hij)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHij() at line" << __LINE__ << std::endl;
#endif

    SquareSubMatrix<double> submat(kbpsi->computeHvnlMatrix(ions));

    SquareSubMatrix2DistMatrix* ss2dm = SquareSubMatrix2DistMatrix::instance();
    ss2dm->accumulate(submat, hij, 0.);

    // add local Hamiltonian part to phi^T*H*phi
    hamiltonian_->addHlocal2matrix(orbitals_i, orbitals_j, hij, false);
}
#endif

template <class OrbitalsType>
void MGmol<OrbitalsType>::computeHij(OrbitalsType& orbitals_i,
    OrbitalsType& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi,
    ProjectedMatricesInterface* projmatrices)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHij() at line " << __LINE__ << std::endl;
#endif

    kbpsi->computeHvnlMatrix(ions, projmatrices);

    // add local Hamiltonian part to phi^T*H*phi
    hamiltonian_->addHlocalij(orbitals_i, orbitals_j, projmatrices);
}

template <class OrbitalsType>
void MGmol<OrbitalsType>::getKBPsiAndHij(OrbitalsType& orbitals_i,
    OrbitalsType& orbitals_j, Ions& ions, KBPsiMatrixSparse* kbpsi,
    ProjectedMatricesInterface* projmatrices)
{
    kbpsi->computeAll(ions, orbitals_i);

    projmatrices->clearSparseH();

    computeHij(orbitals_i, orbitals_j, ions, kbpsi, projmatrices);

    projmatrices->setHiterativeIndex(orbitals_j.getIterativeIndex(),
        hamiltonian_->potential().getIterativeIndex());
}

template <class OrbitalsType>
void MGmol<OrbitalsType>::getKBPsiAndHij(OrbitalsType& orbitals, Ions& ions)
{
    getKBPsiAndHij(
        orbitals, orbitals, ions, g_kbpsi_.get(), proj_matrices_.get());
}

template <class OrbitalsType>
void MGmol<OrbitalsType>::computeHnlPhiAndAdd2HPhi(Ions& ions,
    OrbitalsType& phi, OrbitalsType& hphi, const KBPsiMatrixSparse* const kbpsi)
{
    // H_nl
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHnlPhiAndAdd2HPhi()" << std::endl;
#endif

    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    const int numpt        = mygrid.size();

    const std::vector<std::vector<int>>& gid(phi.getOverlappingGids());
    const short ncolors = gid[0].size();

    {
        using memory_space_type = typename OrbitalsType::memory_space_type;
        // compute Hnl*phi
        if (ct.Mehrstellen())
        {
            pb::GridFuncVector<ORBDTYPE, MemorySpace::Host> gfv(
                mygrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2], gid);
            std::vector<ORBDTYPE> work(numpt * ncolors);
            for (short icolor = 0; icolor < ncolors; icolor++)
            {
                get_vnlpsi(ions, gid, icolor, kbpsi, work.data());
                gfv.assign(icolor, work.data());
            }

            gfv.rhs_4th_Mehr1(work.data());

            // compute B*Hnl*phi and add it to H*phi
            for (int icolor = 0; icolor < ncolors; icolor++)
            {
                // Add the contribution of the non-local potential to H phi
                ORBDTYPE* hpsi           = hphi.getPsi(icolor);
                ORBDTYPE* hpsi_host_view = MemorySpace::Memory<ORBDTYPE,
                    memory_space_type>::allocate_host_view(numpt);
                MemorySpace::Memory<ORBDTYPE,
                    memory_space_type>::copy_view_to_host(hpsi, numpt,
                    hpsi_host_view);

                LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(numpt,
                    (ORBDTYPE)1., work.data() + numpt * icolor, hpsi_host_view);

                MemorySpace::Memory<ORBDTYPE,
                    memory_space_type>::copy_view_to_dev(hpsi_host_view, numpt,
                    hpsi);
                MemorySpace::Memory<ORBDTYPE,
                    memory_space_type>::free_host_view(hpsi_host_view);
            }
        }
        else // no Mehrstellen
        {
            {
                ORBDTYPE* hnl = new ORBDTYPE[numpt];
                for (short icolor = 0; icolor < ncolors; icolor++)
                {
                    get_vnlpsi(ions, gid, icolor, kbpsi, hnl);

                    ORBDTYPE* hpsi           = hphi.getPsi(icolor);
                    ORBDTYPE* hpsi_host_view = MemorySpace::Memory<ORBDTYPE,
                        memory_space_type>::allocate_host_view(numpt);
                    MemorySpace::Memory<ORBDTYPE,
                        memory_space_type>::copy_view_to_host(hpsi, numpt,
                        hpsi_host_view);

                    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
                        numpt, (ORBDTYPE)1., hnl, hpsi_host_view);

                    MemorySpace::Memory<ORBDTYPE,
                        memory_space_type>::copy_view_to_dev(hpsi_host_view,
                        numpt, hpsi);
                    MemorySpace::Memory<ORBDTYPE,
                        memory_space_type>::free_host_view(hpsi_host_view);
                }
                delete[] hnl;
            }
        }
    }

    hphi.setIterativeIndex(phi.getIterativeIndex());
}

template <class OrbitalsType>
void MGmol<OrbitalsType>::getHpsiAndTheta(Ions& ions, OrbitalsType& phi,
    OrbitalsType& hphi, const KBPsiMatrixSparse* const kbpsi)
{
    get_Hpsi_and_Hij_tm_.start();

    const Potentials& pot = hamiltonian_->potential();

    const int phi_it_index = phi.getIterativeIndex();

#ifdef PRINT_OPERATIONS
    os_ << " getHpsiAndTheta" << std::endl;
#endif

    hamiltonian_->applyLocal(phi.chromatic_number(), phi, hphi);

    // Compute "nstates" columns of matrix
    //  Hij = phi**T * H_loc * phi  and save in sh
    if (proj_matrices_->isHupToDate(
            phi.getIterativeIndex(), pot.getIterativeIndex()))
    {
#ifdef PRINT_OPERATIONS
        if (onpe0)
            os_ << "Hij matrix up to date, no computation necessary"
                << std::endl;
#endif
    }
    else
    {
#ifdef PRINT_OPERATIONS
        if (onpe0) os_ << "build matrix Hij = Phi**T * H * Phi" << std::endl;
#endif

        proj_matrices_->clearSparseH();

        // Compute the contribution of the non-local potential into sh
        kbpsi->computeHvnlMatrix(ions, proj_matrices_.get());

        // add local part of H to sh
        SquareLocalMatrices<MATDTYPE, MemorySpace::Host> slm(
            phi.subdivx(), phi.chromatic_number());

        phi.computeLocalProduct(hphi, slm);
        proj_matrices_->setLocalMatrixElementsHl(slm);

        proj_matrices_->consolidateH();

        energy_->saveVofRho();

        // Compute matrix theta = invB * Hij
        proj_matrices_->updateTheta();

        proj_matrices_->setHiterativeIndex(
            phi_it_index, pot.getIterativeIndex());
    }
    computeHnlPhiAndAdd2HPhi(ions, phi, hphi, kbpsi);
    hphi.setIterativeIndex(phi_it_index);

    get_Hpsi_and_Hij_tm_.stop();
}

template class MGmol<LocGridOrbitals<ORBDTYPE>>;
template class MGmol<ExtendedGridOrbitals<ORBDTYPE>>;
