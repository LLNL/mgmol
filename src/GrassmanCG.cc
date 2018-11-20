// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "MGmol.h"

#include "Control.h"
#include "GrassmanCG.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "KBPsiMatrixSparse.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesInterface.h"
#include "ProjectedMatricesSparse.h"

void GrassmanCG::conjugate()
{
    // compute conjugation
    if (conjugate_)
    {
        // Numerator: compute matG = (MG^T*(G-G_old)) = MG^T*G - MG^T*G_old
        // first compute matG = MG^T*G
        SquareLocalMatrices<MATDTYPE> matG(
            new_grad_->subdivx(), new_grad_->chromatic_number());
        new_pcgrad_->getLocalOverlap(*new_grad_,
            matG); // equivalent to pcgrad.computeLocalProduct(grad, matG);
        // compute trace(S^{-1}*matG
        ProjectedMatrices* projmatrices
            = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
        double alpha = projmatrices->computeTraceInvSmultMat(matG);

        // compute matG = MG^T*G_old
        matG.reset();
        new_pcgrad_->getLocalOverlap(*grad_, matG);
        // subtract trace from alpha
        alpha -= projmatrices->computeTraceInvSmultMat(matG);

        // Denominator: compute matG = ((G_old)^T*MG_old)
        matG.reset();
        grad_->getLocalOverlap(*pcgrad_,
            matG); // equivalent to grad_->computeLocalProduct(*pcgrad_, matG);
        // compute trace(S^{-1}*matG
        alpha /= projmatrices->computeTraceInvSmultMat(matG);

        // compute conjugate direction
        double tau       = max(0., alpha);
        const double one = 1.;
        sdir_->scal(tau);
        sdir_->axpy(one, *new_pcgrad_);

        //       if(onpe0)cout<<"conjugate: alpha = "<<alpha<<", tau =
        //       "<<tau<<endl;
    }
    else
    {
        // initialize history data
        grad_   = new LocGridOrbitals("G", *new_grad_, true);
        pcgrad_ = new LocGridOrbitals("P", *new_pcgrad_, true);
        sdir_   = new LocGridOrbitals("S", *new_pcgrad_, true);
    }
}

double GrassmanCG::computeStepSize(LocGridOrbitals& orbitals)
{
    Control& ct   = *(Control::instance());
    const int dim = ct.numst;

    dist_matrix::DistMatrix<DISTMATDTYPE> work_matrix("work_matrix", dim, dim);

    // Done with conjugation. Now compute direction for phi correction.
    // Operations are done assuming an implicit orthogonalization of the
    // search direction with phi (Z = Zo - phi*S^{-1}*phi^T*Zo).
    //
    // The following operations are needed to compute the step size:
    // lambda = Tr[S^{-1}*Z^T*(-G)] / Tr[S^{-1}(Z^T*H*Z) -
    // S^{-1}Z^T*Z*S^{-1}*Phi^T*H*Phi] Where G = natural gradient.

    // Basic ingredients: For seach direction Zo and function Phi
    // Premultiplication by S^{-1} is done to facilitate reusing the matrix
    // 1. dot product Phi^T*H*Zo
    // 2. dot product S^{-1}*Zo^T*H*Phi
    // 3. dot product Zo^T*H*Zo
    // 4. dot product S^{-1}*Zo^T*Phi
    // 5. dot product S^{-1}*Phi^T*Zo
    // 6. dot product Zo^T*Zo
    // We only need to compute and save partial contributions of these
    // matrices on the local subdomains, and consolidate them when needed
    ProjectedMatrices* projmatrices
        = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
    assert(projmatrices);

    // Compute Phi^T*H*Zo and S^{-1}*Zo^T*H*Phi
    dist_matrix::DistMatrix<DISTMATDTYPE> phiHzMat("phiHzMat", dim, dim);
    computeOrbitalsProdWithH(orbitals, *sdir_, phiHzMat);
    // Compute S^{-1}*Zo^T*H*Phi
    dist_matrix::DistMatrix<DISTMATDTYPE> invSzHphiMat(
        "invSzHphiMat", dim, dim);
    invSzHphiMat.transpose(1.0, phiHzMat, 0.);
    projmatrices->applyInvS(invSzHphiMat);

    // compute Zo^T*H*Zo
    dist_matrix::DistMatrix<DISTMATDTYPE> zHzMat("zHzMat", dim, dim);
    computeOrbitalsProdWithH(*sdir_, zHzMat);

    // Compute S^{_1}*Zo^T*Phi and S^{_1}*Phi^T*Zo
    SquareLocalMatrices<MATDTYPE> ss(
        sdir_->subdivx(), sdir_->chromatic_number());
    sdir_->getLocalOverlap(orbitals, ss);
    dist_matrix::DistMatrix<DISTMATDTYPE> invSzTphiMat(
        "invSzTphiMat", dim, dim);
    ss.fillDistMatrix(invSzTphiMat, sdir_->getOverlappingGids());
    dist_matrix::DistMatrix<DISTMATDTYPE> invSphiTzMat("phiTzMat", dim, dim);
    invSphiTzMat.transpose(1.0, invSzTphiMat, 0.);
    // apply invS
    projmatrices->applyInvS(invSzTphiMat);
    projmatrices->applyInvS(invSphiTzMat);

    // Compute Zo^T*Zo
    ss.reset();
    sdir_->getLocalOverlap(ss);
    dist_matrix::DistMatrix<DISTMATDTYPE> zTzMat("zTzMat", dim, dim);
    ss.fillDistMatrix(zTzMat, sdir_->getOverlappingGids());

    // Now compute Tr[S^{-1}*Z^T*(-G)] = Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*H*Phi]
    // - Tr[S^{-1}*Zo^T*H*Phi] Compute Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*H*Phi];
    double anum = projmatrices->computeTraceMatMultTheta(invSzTphiMat);
    // subtract Tr[S^{-1}*Zo^T*H*Phi];
    anum -= invSzHphiMat.trace();

    // done computing numerator.
    // Now compute denominator: Tr[S^{-1}(Z^T*H*Z) -
    // S^{-1}Z^T*Z*S^{-1}*Phi^T*H*Phi] Compute Tr[S^{-1}(Z^T*H*Z) =
    // Tr[S^{-1}*Zo^T*H*Zo]-2*Tr[S^{-1}*Zo^T*H*Phi*S^{-1}*Phi^T*Zo]+Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*H*Phi*S^{-1}*Phi^T*Zo]
    // First compute Tr[S^{-1}*Zo^T*H*Zo]:
    double denom1 = projmatrices->computeTraceInvSmultMat(zHzMat);
    // subtract 2*Tr[S^{-1}*Zo^T*H*Phi*S^{-1}*Phi^T*Zo]:
    work_matrix.gemm('n', 'n', 1., invSzHphiMat, invSphiTzMat, 0.);
    denom1 -= 2 * work_matrix.trace();
    // add Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*H*Phi*S^{-1}*Phi^T*Zo]:
    dist_matrix::DistMatrix<DISTMATDTYPE> pmat("pmat", dim, dim);
    projmatrices->computeMatMultTheta(invSzTphiMat, pmat);
    work_matrix.gemm('n', 'n', 1., pmat, invSphiTzMat, 0.);
    denom1 += work_matrix.trace();

    // Compute S^{-1}Z^T*Z*S^{-1}*Phi^T*H*Phi] =
    // Tr[S^{-1}*Zo^T*Zo*S^{-1}*Phi^T*H*Phi] -
    // Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*Zo*S^{-1}*Phi^T*H*Phi] first compute
    // Tr[S^{-1}*Zo^T*Zo*S^{-1}*Phi^T*H*Phi]:
    double denom2 = projmatrices->computeTraceInvSmultMatMultTheta(zTzMat);
    // subtract Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*Zo*S^{-1}*Phi^T*H*Phi]:
    pmat.clear();
    projmatrices->computeMatMultTheta(invSphiTzMat, pmat);
    work_matrix.gemm('n', 'n', 1., invSzTphiMat, pmat, 0.);
    denom2 -= work_matrix.trace();

    // compute denominator
    double denom = denom1 - denom2;
    // done computing denominator.
    // compute stepsize lambda
    double lambda = anum / denom;

    //   if(onpe0)cout<<"anum = "<<anum<<"; denom1 = "<<denom1<<"; denom2 =
    //   "<<denom2<<"; denom = "<<denom<<"; lambda = "<<lambda<<endl;

    return lambda;
}

// Compute P^T*H*Q for orbitals1-->P and orbitals2-->Q and return result in mat.
// consolidate flag with either gather data or else return local partial
// contributions
void GrassmanCG::computeOrbitalsProdWithH(LocGridOrbitals& orbitals1,
    LocGridOrbitals& orbitals2, dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
{
    // initialize KBPsiMatrices
    KBPsiMatrixSparse kbpsi_1(hamiltonian_->lapOper());
    kbpsi_1.setup(*ptr2ions_, orbitals1);
    kbpsi_1.computeAll(*ptr2ions_, orbitals1);

    KBPsiMatrixSparse kbpsi_2(hamiltonian_->lapOper());
    kbpsi_2.setup(*ptr2ions_, orbitals2);
    kbpsi_2.computeAll(*ptr2ions_, orbitals2);

    // compute P^T*H*Q (orbitals1=P; orbitals2=Q)
    mgmol_strategy_->computeHij(
        orbitals1, orbitals2, *ptr2ions_, &kbpsi_1, &kbpsi_2, mat);

    return;
}

// Compute P^T*H*P for orbitals1-->P and return result in mat.
// consolidate flag with either gather data or else return local partial
// contributions
void GrassmanCG::computeOrbitalsProdWithH(
    LocGridOrbitals& orbitals, dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
{
    // initialize KBPsiMatrices
    KBPsiMatrixSparse kbpsi(hamiltonian_->lapOper());
    kbpsi.setup(*ptr2ions_, orbitals);
    kbpsi.computeAll(*ptr2ions_, orbitals);

    // compute P^T*H*Q (orbitals1=P; orbitals2=Q)
    mgmol_strategy_->computeHij(orbitals, orbitals, *ptr2ions_, &kbpsi, mat);

    return;
}

// parallel transport of history data
// update G=grad_, MG=pcgrad_ and Zo=sdir_
void GrassmanCG::parallelTransportUpdate(
    const double lambda, LocGridOrbitals& phi)
{
    Control& ct = *(Control::instance());

    //    const double fact = lambda;
    const double fact = 1.;

    // update history data
    LocGridOrbitals* gradptr;
    // update gradient information
    gradptr   = new_grad_;
    new_grad_ = grad_;
    grad_     = gradptr;
    // update preconditioned gradient information
    gradptr     = new_pcgrad_;
    new_pcgrad_ = pcgrad_;
    pcgrad_     = gradptr;
    // compute projection-based parallel transport
    if (ct.parallel_transport)
    {
        // compute G_old = G - lambda*(Phi*S^{-1}*Phi^T*G).*corrmasks
        phi.projectOut(*grad_, fact);
        grad_->applyCorrMask(true);
        // compute MG_old = MG - lambda*(Phi*S^{-1}*Phi^T*MG).*masks
        phi.projectOut(*pcgrad_, fact);
        pcgrad_->applyMask(true);
        // update preconditioned search direction information
        phi.projectOut(*sdir_, fact);
        sdir_->applyMask(true);
    }

    // Reset iterative index for current (transported) search direction and
    // gradient This is especially necessary for the search direction, in order
    // to ensure correct computations with the Hamiltonian and wrt phi.
    sdir_->setIterativeIndex(-10);
    grad_->setIterativeIndex(-10);
    pcgrad_->setIterativeIndex(-10);
#if 0
    // test if projection is now 0
    dist_matrix::DistMatrix<DISTMATDTYPE> tmatrix(phi.product(*sdir_));
    if( onpe0 )
        (*MPIdata::sout)<<"GrassmanCG::projectOut(), Product after projection:"<<endl;
    tmatrix.print((*MPIdata::sout),0,0,5,5);
    double trace = tmatrix.trace();
    if(onpe0)cout<<"Trace of projection test = "<<trace<<endl;
#endif

    return;
}
