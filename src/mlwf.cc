// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Control.h"
#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "MGmol.h"
#include "MLWFTransform.h"
#include "Mesh.h"
#include "NOLMOTransform.h"
#include "ProjectedMatrices.h"
#include "ReplicatedWorkSpace.h"
#include "SinCosOps.h"
#include "blas3_c.h"
#include "mputils.h"

#include <vector>
using namespace std;

Timer get_NOLMO_tm("get_NOLMO");
Timer get_MLWF_tm("get_MLWF");

void dtrsm_c(const char side, const char uplo, const char transa,
    const char diag, const int m, const int n, const double alpha,
    const double* const a, const int lda, double* const b, const int ldb)
{
    DTRSM(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

void distributeColumns(
    vector<DISTMATDTYPE>& vmm, dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
{
    assert(mat.nprow() == 1);

    const int nst   = mat.m();
    const int nb    = mat.nb();
    const int count = mat.nloc() * nst;
    // mycol =0 if mat not distributed (duplicated)
    const int mycol = mat.mycol();
    mat.assign(&vmm[mycol * nb * nst], count);
}

template <class T>
int MGmol<T>::getMLWF(MLWFTransform& mlwft, T& orbitals, T& work_orbitals,
    const double dd, const bool apply_flag)
{
    Control& ct = *(Control::instance());
    assert(!ct.isLocMode());

    const int numst = ct.numst;
    if (numst == 0) return 0;

    get_MLWF_tm.start();

    work_orbitals.copyDataFrom(orbitals);

    // orthonormalize work_orbitals before getting sin and cos matrices
    work_orbitals.orthonormalizeLoewdin(false);

    vector<int> overlap(numst * numst);
    int icount = 0;
    for (int st1 = 0; st1 < numst; st1++)
    {
        overlap[st1 + numst * st1] = 1;
        for (int st2 = 0; st2 < st1; st2++)
        {
            double d                   = lrs_->distanceBetweenCenters(st1, st2);
            overlap[st1 + numst * st2] = overlap[st2 + numst * st1]
                = (int)(d < dd);
            icount += overlap[st2 + numst * st1];
        }
    }
    if (onpe0 && ct.verbose > 1)
        os_ << "getMLWF(): Wannier rotations for " << icount << " pairs"
            << endl;
    mlwft.setia(overlap);

    vector<vector<double>> sincos;
    sincos.resize(6);
    for (int i = 0; i < 6; i++)
        sincos[i].resize(numst * numst);

    SinCosOps<T>::compute(work_orbitals, sincos);

    mlwft.distributeColumnsR(sincos);
    // for(int i=0;i<6;i++)
    //   distributeColumns(sincos[i],mlwft.r(i));
    mlwft.compute_transform(200, 1.e-6);

    // mlwft.printTransformationMatrix();

    if (apply_flag)
    {
        if (onpe0 && ct.verbose > 1) os_ << "Apply Wannier rotations" << endl;

        assert(&orbitals != &work_orbitals);

        // orthogonal transformation of orbitals
        work_orbitals.multiply_by_matrix(&mlwft.mat()[0], orbitals);
    }

    orbitals.incrementIterativeIndex();

    get_MLWF_tm.stop();

    return 0;
}

template <class T>
int MGmol<T>::getMLWF2states(
    const int st1, const int st2, T& orbitals, T& work_orbitals)
{
    get_MLWF_tm.start();

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    Vector3D ll(mygrid.ll(0), mygrid.ll(1), mygrid.ll(2));
    Vector3D origin(mygrid.origin(0), mygrid.origin(1), mygrid.origin(2));
    MLWFTransform mlwft(2, origin, ll);

    Control& ct = *(Control::instance());

    work_orbitals.copyDataFrom(orbitals);

    // orthonormalize work_orbitals to get sin and cos matrices
    work_orbitals.orthonormalize2states(st1, st2);

    vector<int> overlap(4, 1);
    if (onpe0 && ct.verbose > 1)
        os_ << "getMLWF2states(): Wannier rotations for 1 pair" << endl;
    mlwft.setia(overlap);

    vector<vector<double>> sincos;
    sincos.resize(6);
    for (int i = 0; i < 6; i++)
        sincos[i].resize(4);

    SinCosOps<T>::compute2states(work_orbitals, sincos, st1, st2);

    for (int i = 0; i < 6; i++)
        distributeColumns(sincos[i], mlwft.r(i));
    mlwft.compute_transform(200, 1.e-6);

    if (onpe0 && ct.verbose > 1)
    {
        os_ << "Apply Wannier rotations for 2 states" << endl;
        if (ct.verbose > 2) mlwft.printTransformationMatrix();
    }
    assert(&orbitals != &work_orbitals);

    // orthogonal transformation of orbitals
    work_orbitals.multiplyByMatrix2states(st1, st2, &mlwft.mat()[0], orbitals);

    orbitals.incrementIterativeIndex();

    get_MLWF_tm.stop();

    return 0;
}

template <class T>
void MGmol<T>::wftransform(T* orbitals, T* work_orbitals, Ions& ions)
{
    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    assert(!ct.isLocMode());

    bool createMLWF = false;
    bool createNOOT = false;
    double dd       = 1000.;
    if (ct.isLocMode()) dd = 0.1;
    const int numst = ct.numst;
    Vector3D origin(mygrid.origin(0), mygrid.origin(1), mygrid.origin(2));
    Vector3D ll(mygrid.ll(0), mygrid.ll(1), mygrid.ll(2));

    MLWFTransform* mlwt  = nullptr;
    NOLMOTransform* noot = nullptr;
    if (numst > 0)
    {
        mlwt = new MLWFTransform(numst, origin, ll);

        createMLWF          = true;
        bool apply_rotation = (ct.wannier_transform_type >= 2);
        getMLWF(*mlwt, *orbitals, *work_orbitals, dd, apply_rotation);

        // need to recompute matrices, even if transform was not applied
        // because orbitals and work_orbitals share projected matrices
        resetProjectedMatricesAndDM(*orbitals, ions);

        if (ct.verbose > 1 || !ct.AtomsMove())
            mlwt->printCentersAndSpreads(os_);
        if (ct.wannier_transform_type == 3)
        {
            if (onpe0) os_ << "NO orbitals centers and spread" << endl;
            noot       = new NOLMOTransform(numst, origin, ll);
            createNOOT = true;

            get_NOLMO(*noot, *orbitals, *work_orbitals, dd, false);
            if (ct.verbose > 1 || !ct.AtomsMove())
                noot->printCentersAndSpreads(os_);
        }
    }

    // print transformation matrix
    if (ct.getOrbitalsType() == OrbitalsType::Eigenfunctions && !ct.AtomsMove())
    {
        mlwt->printTransform();
    }

    if (createMLWF)
    {
        assert(mlwt != NULL);
        delete mlwt;
    }
    if (createNOOT)
    {
        assert(noot != NULL);
        delete noot;
    }
}

template <class T>
int MGmol<T>::get_NOLMO(NOLMOTransform& noot, T& orbitals, T& work_orbitals,
    const double dd, const bool apply_flag)
{
    get_NOLMO_tm.start();

    Control& ct = *(Control::instance());

    work_orbitals.copyDataFrom(orbitals);

    int numst = ct.numst;
    vector<int> overlap(numst * numst);
    int icount = 0;
    for (int st1 = 0; st1 < numst; st1++)
    {
        overlap[st1 + numst * st1] = 1;
        for (int st2 = 0; st2 < st1; st2++)
        {
            double d                   = lrs_->distanceBetweenCenters(st1, st2);
            overlap[st1 + numst * st2] = overlap[st2 + numst * st1]
                = (int)(d < dd);
            icount += overlap[st2 + numst * st1];
        }
    }
    if (onpe0 && ct.verbose > 1)
        os_ << "get_NOLMO(): for " << icount << " pairs" << endl;
    noot.setia(overlap);

    // compute sin/cos matrices one dimension at a time
    vector<vector<double>> sincos;
    sincos.resize(2);
    for (short i = 0; i < 2; i++)
        sincos[i].resize(numst * numst);

    for (short d = 0; d < 3; d++)
    {
        SinCosOps<T>::compute1D(work_orbitals, sincos, d);
        for (short i = 0; i < 2; i++)
            distributeColumns(sincos[i], noot.r(i + 2 * d));

        SinCosOps<T>::computeSquare1D(work_orbitals, sincos, d);
        for (short i = 0; i < 2; i++)
            distributeColumns(sincos[i], noot.b(i + 2 * d));
    }
    sincos.clear();

    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    assert(projmatrices);

    noot.init_transform(projmatrices->getLS(), true);
    noot.compute_transform(100, 1.e-6);

    if (apply_flag)
    {
        if (onpe0) os_ << "Apply NOLMO rotations" << endl;

        assert(&orbitals != &work_orbitals);

        DISTMATDTYPE* a = &noot.mat()[0];

        // include initial a0 into transformation a
        ReplicatedWorkSpace<DISTMATDTYPE>& wspace(
            ReplicatedWorkSpace<DISTMATDTYPE>::instance());
        wspace.initSquareMatrix(projmatrices->getLS());
        wspace.setUpperTriangularSquareMatrixToZero();

        dtrsm_c('l', 'l', 'n', 'n', numst, numst, 1., wspace.square_matrix(),
            numst, a, numst);

        // normalize columns of a to get normalized functions
        vector<double> alphan(numst);
        for (int p = 0; p < numst; p++)
        {
            alphan[p] = 0.;
            for (int q = 0; q < numst; q++)
                alphan[p] += (a[q + numst * p] * a[q + numst * p]);
            // os_ << " norm:   " << sqrt(alphan[p]) << endl;
            alphan[p] = sqrt(alphan[p]);
            for (int q = 0; q < numst; q++)
                a[q + numst * p] /= alphan[p];
        }

        // nonorthogonal transformation of orbitals
        // result in orbitals
        work_orbitals.multiply_by_matrix(a, orbitals);
        orbitals.incrementIterativeIndex();

        dist_matrix::DistMatrix<DISTMATDTYPE> u_dis("Udis", numst, numst);
        u_dis.init(a, numst);

        ProjectedMatrices* proj_matrices
            = dynamic_cast<ProjectedMatrices*>(orbitals.proj_matrices());
        proj_matrices->rotateAll(u_dis, false);
        // orbitals.rotateSubMatrices(u_dis);
    }

    get_NOLMO_tm.stop();

    return 0;
}

template class MGmol<LocGridOrbitals>;
template class MGmol<ExtendedGridOrbitals>;
