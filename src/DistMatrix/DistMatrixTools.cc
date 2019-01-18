// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DistMatrixTools.h"
#include "MPIdata.h"
#include "MatricesBlacsContext.h"

#include <vector>
using namespace std;

void getProcrustesTransform(dist_matrix::DistMatrix<DISTMATDTYPE>& q,
    dist_matrix::DistMatrix<DISTMATDTYPE>& yyt)
{
    if (onpe0) (*MPIdata::sout) << "getProcrustesTransform()" << endl;
    MatricesBlacsContext& mbc(MatricesBlacsContext::instance());
    const dist_matrix::BlacsContext& bc = *mbc.bcxt();
    const int nst                       = q.m();
    const DISTMATDTYPE tol              = 0.5;

    dist_matrix::DistMatrix<DISTMATDTYPE> z("z", bc, nst, nst);
    dist_matrix::DistMatrix<DISTMATDTYPE> yt("yt", bc, nst, nst);

    vector<DISTMATDTYPE> eigenvalues(nst);

    // get z and yt
    q.gesvd('V', 'V', eigenvalues, z, yt);
    if (onpe0)
    {
        (*MPIdata::sout) << "Procrustes: Small singular values: ";
        for (int i = 0; i < nst; i++)
            if (eigenvalues[i] < tol)
                (*MPIdata::sout) << "  " << eigenvalues[i];
        (*MPIdata::sout) << endl;
    }

    // get q=z*yt
    q.gemm('n', 'n', 1., z, yt, 0.);

    // now calculate yyt matrix
    vector<DISTMATDTYPE> diag(nst);
    for (int i = 0; i < nst; i++)
        diag[i] = (eigenvalues[i] < tol) ? 0. : 1.;
    dist_matrix::DistMatrix<DISTMATDTYPE> gamma(
        "Gamma", bc, &diag[0], nst, nst);
    z.gemm('t', 'n', 1., yt, gamma, 0.);
    yyt.gemm('n', 'n', 1., z, yt, 0.);
}

// Find
// rotation = arg min_U ||psi*U-phi||_F
// by computing U=(O*O^T)^{-1/2}*O
void getAlignFrobeniusTransform(dist_matrix::DistMatrix<DISTMATDTYPE>& omatrix,
    dist_matrix::DistMatrix<DISTMATDTYPE>& rotation)
{
    MatricesBlacsContext& mbc(MatricesBlacsContext::instance());
    const dist_matrix::BlacsContext& bc = *mbc.bcxt();
    const int nst                       = omatrix.m();

    // get z=o*o^t
    dist_matrix::DistMatrix<DISTMATDTYPE> zmat("z", bc, nst, nst);
    zmat.gemm('n', 't', 1., omatrix, omatrix, 0.);

    dist_matrix::DistMatrix<DISTMATDTYPE> vmat("v", bc, nst, nst);
    vector<DISTMATDTYPE> eigenvalues(nst);

    zmat.syev('V', 'u', eigenvalues, vmat);

    for (int i = 0; i < nst; i++)
        eigenvalues[i] = 1. / sqrt(eigenvalues[i]);

    dist_matrix::DistMatrix<DISTMATDTYPE> gamma(
        "Gamma", bc, &eigenvalues[0], nst, nst);

    dist_matrix::DistMatrix<DISTMATDTYPE> work("w", bc, nst, nst);
    // work = vmat*gamma with gamma symmetric
    work.symm('r', 'l', 1., gamma, vmat, 0.);

    // zmat = work * vmat^T
    zmat.gemm('n', 't', 1., work, vmat, 0.);

    rotation.gemm('n', 'n', 1., zmat, omatrix, 0.);
}
