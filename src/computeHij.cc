// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cassert>

#include "Control.h"
#include "DistMatrix.h"
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
#include "SquareSubMatrix2DistMatrix.h"

template <>
void MGmol<LocGridOrbitals>::addHlocal2matrix(LocGridOrbitals& orbitalsi,
    LocGridOrbitals& orbitalsj, VariableSizeMatrix<sparserow>& mat)
{
    computeHij_tm_.start();

#if DEBUG
    os_ << " addHlocal2matrix()" << endl;
#endif

    hamiltonian_->addHlocal2matrix(orbitalsi, orbitalsj, mat, true);

    computeHij_tm_.stop();
}

template <>
void MGmol<LocGridOrbitals>::addHlocal2matrix(LocGridOrbitals& orbitalsi,
    LocGridOrbitals& orbitalsj, dist_matrix::DistMatrix<double>& H)
{
    computeHij_tm_.start();

#if DEBUG
    os_ << " addHlocal2matrix()" << endl;
#endif

    hamiltonian_->addHlocal2matrix(orbitalsi, orbitalsj, H);

    computeHij_tm_.stop();
}

template <>
template <>
void MGmol<LocGridOrbitals>::computeHij(LocGridOrbitals& orbitals_i,
    LocGridOrbitals& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi_i,
    const KBPsiMatrixSparse* const kbpsi_j, VariableSizeMatrix<sparserow>& mat,
    const bool consolidate)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHij()" << endl;
#endif

    // compute phi_i^T*Hnl*Phi_j
    kbpsi_i->computeHvnlMatrix(kbpsi_j, ions, mat);

    // add local Hamiltonian part to phi_i^T*H*phi_j
    addHlocal2matrix(orbitals_i, orbitals_j, mat);

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
void MGmol<LocGridOrbitals>::computeHij(LocGridOrbitals& orbitals_i,
    LocGridOrbitals& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi, VariableSizeMatrix<sparserow>& mat,
    const bool consolidate)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHij()" << endl;
#endif

    kbpsi->computeHvnlMatrix(ions, mat);

    // add local Hamiltonian part to phi^T*H*phi
    addHlocal2matrix(orbitals_i, orbitals_j, mat);

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

template <class T>
void MGmol<T>::computeHij_private(T& orbitals_i, T& orbitals_j,
    const Ions& ions, const KBPsiMatrixSparse* const kbpsi_i,
    const KBPsiMatrixSparse* const kbpsi_j,
    dist_matrix::DistMatrix<DISTMATDTYPE>& hij)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHij()" << endl;
#endif

    hij.clear();

    SquareSubMatrix<double> submat(kbpsi_i->computeHvnlMatrix(kbpsi_j, ions));

    SquareSubMatrix2DistMatrix* ss2dm = SquareSubMatrix2DistMatrix::instance();
    ss2dm->accumulate(submat, hij, 0.);

    // add local Hamiltonian part to phi^T*H*phi
    addHlocal2matrix(orbitals_i, orbitals_j, hij);
}

template <>
template <>
void MGmol<LocGridOrbitals>::computeHij(LocGridOrbitals& orbitals_i,
    LocGridOrbitals& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi,
    const KBPsiMatrixSparse* const kbpsi_j,
    dist_matrix::DistMatrix<DISTMATDTYPE>& hij, const bool consolidate)
{
    (void)consolidate;

    computeHij_private(orbitals_i, orbitals_j, ions, kbpsi, kbpsi_j, hij);
}

template <>
template <>
void MGmol<ExtendedGridOrbitals>::computeHij(ExtendedGridOrbitals& orbitals_i,
    ExtendedGridOrbitals& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi,
    const KBPsiMatrixSparse* const kbpsi_j,
    dist_matrix::DistMatrix<DISTMATDTYPE>& hij, const bool consolidate)
{
    (void)consolidate;

    computeHij_private(orbitals_i, orbitals_j, ions, kbpsi, kbpsi_j, hij);
}

template <class T>
void MGmol<T>::computeHij(T& orbitals_i, T& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi, dist_matrix::DistMatrix<double>& hij,
    const bool consolidate)
{
    (void)consolidate;

    computeHij_private(orbitals_i, orbitals_j, ions, kbpsi, hij);
}

template <class T>
void MGmol<T>::computeHij_private(T& orbitals_i, T& orbitals_j,
    const Ions& ions, const KBPsiMatrixSparse* const kbpsi,
    dist_matrix::DistMatrix<double>& hij)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHij()" << endl;
#endif

    SquareSubMatrix<double> submat(kbpsi->computeHvnlMatrix(ions));

    SquareSubMatrix2DistMatrix* ss2dm = SquareSubMatrix2DistMatrix::instance();
    ss2dm->accumulate(submat, hij, 0.);

    // add local Hamiltonian part to phi^T*H*phi
    addHlocal2matrix(orbitals_i, orbitals_j, hij);
}

template <class T>
void MGmol<T>::computeHij(T& orbitals_i, T& orbitals_j, const Ions& ions,
    const KBPsiMatrixSparse* const kbpsi,
    ProjectedMatricesInterface* projmatrices)
{
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHij()" << endl;
#endif

    kbpsi->computeHvnlMatrix(ions, projmatrices);

    // add local Hamiltonian part to phi^T*H*phi
    addHlocalij(orbitals_i, orbitals_j, projmatrices);

    projmatrices->consolidateH();
}

template <class T>
void MGmol<T>::getKBPsiAndHij(T& orbitals_i, T& orbitals_j, Ions& ions,
    KBPsiMatrixSparse* kbpsi, ProjectedMatricesInterface* projmatrices,
    dist_matrix::DistMatrix<DISTMATDTYPE>& hij)
{
    kbpsi->computeAll(ions, orbitals_i);

    computeHij(orbitals_i, orbitals_j, ions, kbpsi, hij, true);

    projmatrices->setHiterativeIndex(orbitals_j.getIterativeIndex(),
        hamiltonian_->potential().getIterativeIndex());
}

template <class T>
void MGmol<T>::getKBPsiAndHij(T& orbitals_i, T& orbitals_j, Ions& ions,
    KBPsiMatrixSparse* kbpsi, ProjectedMatricesInterface* projmatrices)
{
    kbpsi->computeAll(ions, orbitals_i);

    projmatrices->clearSparseH();

    computeHij(orbitals_i, orbitals_j, ions, kbpsi, projmatrices);

    projmatrices->setHiterativeIndex(orbitals_j.getIterativeIndex(),
        hamiltonian_->potential().getIterativeIndex());
}

template <class T>
void MGmol<T>::getKBPsiAndHij(T& orbitals, Ions& ions, KBPsiMatrixSparse* kbpsi,
    dist_matrix::DistMatrix<DISTMATDTYPE>& hij)
{
    getKBPsiAndHij(orbitals, orbitals, ions, kbpsi, proj_matrices_, hij);
}

template <class T>
void MGmol<T>::getKBPsiAndHij(T& orbitals, Ions& ions)
{
    getKBPsiAndHij(orbitals, orbitals, ions, g_kbpsi_, proj_matrices_);
}

template <class T>
void MGmol<T>::computeHnlPhiAndAdd2HPhi(
    Ions& ions, T& phi, T& hphi, const KBPsiMatrixSparse* const kbpsi)
{
    // H_nl
#ifdef PRINT_OPERATIONS
    if (onpe0) os_ << "computeHnlPhiAndAdd2HPhi()" << endl;
#endif

    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    const int numpt        = mygrid.size();

    const std::vector<std::vector<int>>& gid(phi.getOverlappingGids());
    const short ncolors = gid[0].size();

    {
        // compute Hnl*phi
        if (ct.Mehrstellen())
        {
            pb::GridFuncVector<ORBDTYPE> gfv(
                true, mygrid, ct.bc[0], ct.bc[1], ct.bc[2], gid);
            std::vector<ORBDTYPE> work(numpt);
            for (short icolor = 0; icolor < ncolors; icolor++)
            {
                get_vnlpsi(ions, gid, icolor, kbpsi, work.data());
                gfv.assign(icolor, work.data());
            }
            gfv.trade_boundaries();

            // compute B*Hnl*phi and add it to H*phi
            for (int icolor = 0; icolor < ncolors; icolor++)
            {
                hamiltonian_->lapOper()->rhs(gfv.func(icolor), work.data());

                // Add the contribution of the non-local potential to H phi
                ORBDTYPE* hpsi = hphi.getPsi(icolor);
                MPaxpy(numpt, 1., work.data(), hpsi);
            }
        }
        else // no Mehrstellen
        {
#pragma omp parallel
            {
                ORBDTYPE* hnl = new ORBDTYPE[numpt];
#pragma omp for
                for (short icolor = 0; icolor < ncolors; icolor++)
                {
                    get_vnlpsi(ions, gid, icolor, kbpsi, hnl);
                    MPaxpy(numpt, 1., hnl, hphi.getPsi(icolor));
                }

                delete[] hnl;
            }
        }
    }

    hphi.setIterativeIndex(phi.getIterativeIndex());
}

template <class T>
void MGmol<T>::addHlocal2matrix(
    T& orbitalsi, T& orbitalsj, dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
{
    computeHij_tm_.start();

#if DEBUG
    os_ << " addHlocal2matrix()" << endl;
#endif

    // add local H to DistMatrix
    hamiltonian_->addHlocal2matrix(orbitalsi, orbitalsj, mat);

    computeHij_tm_.stop();
}

template <class T>
void MGmol<T>::addHlocalij(
    T& orbitalsi, T& orbitalsj, ProjectedMatricesInterface* pmat)
{
    computeHij_tm_.start();

#if DEBUG
    os_ << " addHlocalij()" << endl;
#endif

    hamiltonian_->addHlocalij(orbitalsi, orbitalsj, pmat);

    computeHij_tm_.stop();
}

template <class T>
void MGmol<T>::getHpsiAndTheta(Ions& ions, T& phi, T& hphi)
{
    getHpsiAndTheta(ions, phi, hphi, g_kbpsi_);
}

template <class T>
void MGmol<T>::getHpsiAndTheta(
    Ions& ions, T& phi, T& hphi, const KBPsiMatrixSparse* const kbpsi)
{
    get_Hpsi_and_Hij_tm_.start();

    const Potentials& pot = hamiltonian_->potential();

    const int phi_it_index = phi.getIterativeIndex();

#if DEBUG
    os_ << " getHpsiAndTheta" << endl;
#endif

    hphi.assign(hamiltonian_->applyLocal(phi));

    // Compute "nstates" columns of matrix
    //  Hij = phi**T * H_loc * phi  and save in sh
    if (proj_matrices_->isHupToDate(
            phi.getIterativeIndex(), pot.getIterativeIndex()))
    {
#ifdef PRINT_OPERATIONS
        if (onpe0)
            os_ << "Hij matrix up to date, no computation necessary" << endl;
#endif
    }
    else
    {
#ifdef PRINT_OPERATIONS
        if (onpe0) os_ << "build matrix Hij = Phi**T * H * Phi" << endl;
#endif

        proj_matrices_->clearSparseH();

        // Compute the contribution of the non-local potential into sh
        kbpsi->computeHvnlMatrix(ions, proj_matrices_);

        // add local part of H to sh
        SquareLocalMatrices<MATDTYPE> slh(
            phi.subdivx(), phi.chromatic_number());

        phi.computeLocalProduct(hphi, slh);

        proj_matrices_->setLocalMatrixElementsHl(slh);

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

template class MGmol<LocGridOrbitals>;
template class MGmol<ExtendedGridOrbitals>;
