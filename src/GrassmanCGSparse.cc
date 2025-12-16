// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "MGmol.h"

#include "Control.h"
#include "GrassmanCGSparse.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "KBPsiMatrixSparse.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesInterface.h"
#include "ProjectedMatricesSparse.h"

template <class T>
void GrassmanCGSparse<T>::conjugate()
{
    // compute conjugation
    // We tranpose the dot product operations where necessary. This is done
    // for efficiency, since matrix operations which use VariableSizeMatrix
    // format are done on the rows of the matrix (or assumes symmetry of
    // matrix). If the dot products are transposed ahead of time, then there is
    // no need to transpose the resulting matrices prior to performing
    // operations using the VariableSizeMatrix format.
    if (GrassmanLineMinimization<T>::conjugate_)
    {
        ProjectedMatricesSparse* projmatrices
            = dynamic_cast<ProjectedMatricesSparse*>(
                GrassmanLineMinimization<T>::proj_matrices_);

        // compute G-G_old
        const double m_one = -1.;
        GrassmanLineMinimization<T>::new_grad_->axpy(
            m_one, *GrassmanLineMinimization<T>::grad_);
        // Numerator: compute matG = (MG^T*(G-G_old))^T = G^T*MG - G_old^T*MG
        // first compute matG = G^T*MG
        SquareLocalMatrices<MATDTYPE, MemorySpace::Host> matG(
            GrassmanLineMinimization<T>::new_grad_->subdivx(),
            GrassmanLineMinimization<T>::new_grad_->chromatic_number());
        GrassmanLineMinimization<T>::new_grad_->getLocalOverlap(
            *GrassmanLineMinimization<T>::new_pcgrad_, matG);
        // compute trace(S^{-1}*matG
        double alpha = projmatrices->computeTraceInvSmultMat(matG, true);

        // compute matG = G_old^T*MG
        matG.reset();
        GrassmanLineMinimization<T>::grad_->getLocalOverlap(
            *GrassmanLineMinimization<T>::new_pcgrad_, matG);
        // subtract trace from alpha
        alpha -= projmatrices->computeTraceInvSmultMat(matG, true);

        // Denominator: compute matG = ((G_old)^T*MG_old)^T = MG_old^T*G_old
        matG.reset();
        GrassmanLineMinimization<T>::pcgrad_->getLocalOverlap(
            *GrassmanLineMinimization<T>::grad_, matG);
        // compute trace(S^{-1}*matG
        alpha /= projmatrices->computeTraceInvSmultMat(matG, true);

        // compute conjugate direction
        double tau       = std::max(0., alpha);
        const double one = 1.;
        GrassmanLineMinimization<T>::sdir_->scal(tau);
        GrassmanLineMinimization<T>::sdir_->axpy(
            one, *GrassmanLineMinimization<T>::new_pcgrad_);
    }
    else
    {
        // initialize history data
        GrassmanLineMinimization<T>::grad_
            = new T("G", *GrassmanLineMinimization<T>::new_grad_, true);
        GrassmanLineMinimization<T>::pcgrad_
            = new T("P", *GrassmanLineMinimization<T>::new_pcgrad_, true);
        GrassmanLineMinimization<T>::sdir_
            = new T("S", *GrassmanLineMinimization<T>::new_pcgrad_, true);
    }
}

template <class T>
double GrassmanCGSparse<T>::computeStepSize(T& orbitals)
{

    // Done with conjugation. Now compute direction for phi correction.
    // Operations are done assuming an implicit orthogonalization of the
    // search direction with phi (Z = Zo - phi*S^{-1}*phi^T*Zo).
    //
    // The following operations are needed to compute the step size:
    // lambda = Tr[S^{-1}*Z^T*G] / Tr[S^{-1}(Z^T*H*Z) -
    // S^{-1}Z^T*Z*S^{-1}*Phi^T*H*Phi] Where G = natural gradient.

    // Basic ingredients: For seach direction Zo and function Phi
    // 1. dot product Zo^T*H*Phi
    // 2. dot product Phi^T*H*Zo
    // 3. dot product Zo^T*H*Zo
    // 4. dot product Zo^T*Phi
    // 5. dot product Phi^T*Zo
    // 6. dot product Zo^T*Zo
    // We only need to compute and save partial contributions of these
    // matrices on the local subdomains, and consolidate them when needed
    Control& ct = *(Control::instance());

    // Compute Z^T*H*Phi
    VariableSizeMatrix<sparserow> zHphiMat(
        "zHphi", GrassmanLineMinimization<T>::sdir_->chromatic_number());
    computeOrbitalsProdWithH(
        *GrassmanLineMinimization<T>::sdir_, orbitals, zHphiMat, false);
    // Compute Phi^T*H*Z
    VariableSizeMatrix<sparserow> phiHzMat(
        "phiHz", orbitals.chromatic_number());
    computeOrbitalsProdWithH(
        orbitals, *GrassmanLineMinimization<T>::sdir_, phiHzMat, false);
    // compute Zo^T*H*Zo
    VariableSizeMatrix<sparserow> zHzMat(
        "zHz", GrassmanLineMinimization<T>::sdir_->chromatic_number());
    computeOrbitalsProdWithH(
        *GrassmanLineMinimization<T>::sdir_, zHzMat, false);
    // Compute Zo^T*Phi
    SquareLocalMatrices<MATDTYPE, MemorySpace::Host> ss(
        GrassmanLineMinimization<T>::sdir_->subdivx(),
        GrassmanLineMinimization<T>::sdir_->chromatic_number());
    GrassmanLineMinimization<T>::sdir_->getLocalOverlap(orbitals, ss);
    VariableSizeMatrix<sparserow> zTphiMat(
        "zTphi", GrassmanLineMinimization<T>::sdir_->chromatic_number());
    zTphiMat.insertMatrixElements(
        ss, GrassmanLineMinimization<T>::sdir_->getOverlappingGids(), ct.numst);
    // Compute Phi^T*Zo
    ss.reset();
    orbitals.getLocalOverlap(*GrassmanLineMinimization<T>::sdir_, ss);
    VariableSizeMatrix<sparserow> phiTzMat(
        "phiTz", orbitals.chromatic_number());
    zTphiMat.insertMatrixElements(ss, orbitals.getOverlappingGids(), ct.numst);
    // Compute Zo^T*Zo
    ss.reset();
    GrassmanLineMinimization<T>::sdir_->getLocalOverlap(ss);
    VariableSizeMatrix<sparserow> zTzMat(
        "zTz", GrassmanLineMinimization<T>::sdir_->chromatic_number());
    zTzMat.insertMatrixElements(
        ss, GrassmanLineMinimization<T>::sdir_->getOverlappingGids(), ct.numst);

    // Now compute Tr[S^{-1}*Z^T*G] =
    // Tr[S^{-1}*Zo^T*H*Phi]-Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*H*Phi]
    // Tr[S^{-1}*Zo^T*H*Phi] = Tr[S^{-1}*(Phi^T*H*Zo)^T] <-- matrices are stored
    // and accessed row-wise should we use projmatrices_ ... yes. Sq. local
    // matrices cover only local subdomain and may miss some terms (rows)

    // Now compute Tr[S^{-1}*Z^T*(-G)] = Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*H*Phi]
    // - Tr[S^{-1}*Zo^T*H*Phi] Compute Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*H*Phi] =
    //   double anum = proj_matrices_->computeTraceMatMultTheta(invSzTphiMat);
    // subtract Tr[S^{-1}*Zo^T*H*Phi];
    //   anum -= invSzHphiMat.trace();

    double lambda = 0.;

    return lambda;
}

/*
double GrassmanCG::computeStepSize(T& orbitals)
{
    Control& ct = *(Control::instance());
    const int dim = ct.numst;

   dist_matrix::DistMatrix<DISTMATDTYPE> work_matrix("work_matrix", dim, dim);

   // Done with conjugation. Now compute direction for phi correction.
   // Operations are done assuming an implicit orthogonalization of the
   // search direction with phi (Z = Zo - phi*S^{-1}*phi^T*Zo).
   //
   // The following operations are needed to compute the step size:
   // lambda = Tr[S^{-1}*Z^T*(-G)] / Tr[S^{-1}(Z^T*H*Z) -
S^{-1}Z^T*Z*S^{-1}*Phi^T*H*Phi]
   // Where G = natural gradient.

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

    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();

   // Compute Phi^T*H*Zo and S^{-1}*Zo^T*H*Phi
   dist_matrix::DistMatrix<DISTMATDTYPE> phiHzMat("phiHzMat", dim, dim);
   computeOrbitalsProdWithH(orbitals, *sdir_, phiHzMat);
   // Compute S^{-1}*Zo^T*H*Phi
   dist_matrix::DistMatrix<DISTMATDTYPE> invSzHphiMat("invSzHphiMat", dim, dim);
   invSzHphiMat.transpose(1.0, phiHzMat, 0.);
   proj_matrices_->applyInvS(invSzHphiMat);

   // compute Zo^T*H*Zo
   dist_matrix::DistMatrix<DISTMATDTYPE> zHzMat("zHzMat", dim, dim);
   computeOrbitalsProdWithH(*sdir_, zHzMat);

   // Compute S^{_1}*Zo^T*Phi and S^{_1}*Phi^T*Zo
   SquareLocalMatrices<MATDTYPE,MemorySpace::Host> ss(sdir_->subdivx(),
sdir_->chromatic_number()); sdir_->getLocalOverlap(orbitals, ss);
   dist_matrix::DistMatrix<DISTMATDTYPE> invSzTphiMat("invSzTphiMat", dim, dim);
   ss.fillDistMatrix(invSzTphiMat, sdir_->getOverlappingGids());
   dist_matrix::DistMatrix<DISTMATDTYPE> invSphiTzMat("phiTzMat", dim, dim);
   invSphiTzMat.transpose(1.0, invSzTphiMat, 0.);
   // apply invS
   proj_matrices_->applyInvS(invSzTphiMat);
   proj_matrices_->applyInvS(invSphiTzMat);

   // Compute Zo^T*Zo
   ss.reset();
   sdir_->getLocalOverlap(ss);
   dist_matrix::DistMatrix<DISTMATDTYPE> zTzMat("zTzMat", dim, dim);
   ss.fillDistMatrix(zTzMat, sdir_->getOverlappingGids());

   // Now compute Tr[S^{-1}*Z^T*(-G)] = Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*H*Phi] -
Tr[S^{-1}*Zo^T*H*Phi]
   // Compute Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*H*Phi];
   double anum = proj_matrices_->computeTraceMatMultTheta(invSzTphiMat);
   //subtract Tr[S^{-1}*Zo^T*H*Phi];
   anum -= invSzHphiMat.trace();

   //done computing numerator.
   //Now compute denominator: Tr[S^{-1}(Z^T*H*Z) -
S^{-1}Z^T*Z*S^{-1}*Phi^T*H*Phi]
   //Compute Tr[S^{-1}(Z^T*H*Z) =
Tr[S^{-1}*Zo^T*H*Zo]-2*Tr[S^{-1}*Zo^T*H*Phi*S^{-1}*Phi^T*Zo]+Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*H*Phi*S^{-1}*Phi^T*Zo]
   //First compute Tr[S^{-1}*Zo^T*H*Zo]:
   double denom1 = proj_matrices_->computeTraceInvSmultMat(zHzMat);
   //subtract 2*Tr[S^{-1}*Zo^T*H*Phi*S^{-1}*Phi^T*Zo]:
   work_matrix.gemm('n','n',1.,invSzHphiMat,invSphiTzMat,0.);
   denom1 -=2*work_matrix.trace();
   //add Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*H*Phi*S^{-1}*Phi^T*Zo]:
   dist_matrix::DistMatrix<DISTMATDTYPE> pmat("pmat", dim, dim);
   proj_matrices_->computeMatMultTheta(invSzTphiMat,pmat);
   work_matrix.gemm('n','n',1.,pmat,invSphiTzMat,0.);
   denom1 += work_matrix.trace();

   //Compute S^{-1}Z^T*Z*S^{-1}*Phi^T*H*Phi] =
Tr[S^{-1}*Zo^T*Zo*S^{-1}*Phi^T*H*Phi] -
Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*Zo*S^{-1}*Phi^T*H*Phi]
   //first compute Tr[S^{-1}*Zo^T*Zo*S^{-1}*Phi^T*H*Phi]:
   double denom2 = proj_matrices_->computeTraceInvSmultMatMultTheta(zTzMat);
   //subtract Tr[S^{-1}*Zo^T*Phi*S^{-1}*Phi^T*Zo*S^{-1}*Phi^T*H*Phi]:
   pmat.clear();
   proj_matrices_->computeMatMultTheta(invSphiTzMat,pmat);
   work_matrix.gemm('n','n',1.,invSzTphiMat,pmat,0.);
   denom2 -= work_matrix.trace();

   // compute denominator
   double denom = denom1 - denom2;
   //done computing denominator.
   // compute stepsize lambda
   double lambda= anum/denom;

//   if(onpe0)cout<<"anum = "<<anum<<"; denom1 = "<<denom1<<"; denom2 =
"<<denom2<<"; denom = "<<denom<<"; lambda = "<<lambda<<endl;

   return lambda;
}


*/

/////////////////////////////////////

// Compute P^T*H*Q for orbitals1-->P and orbitals2-->Q and return result in mat.
// consolidate flag with either gather data or else return local partial
// contributions
template <class T>
void GrassmanCGSparse<T>::computeOrbitalsProdWithH(T& orbitals1, T& orbitals2,
    VariableSizeMatrix<sparserow>& mat, const bool consolidate)
{
    // initialize KBPsiMatrices
    KBPsiMatrixSparse kbpsi_1(
        GrassmanLineMinimization<T>::hamiltonian_->lapOper());
    kbpsi_1.setup(*GrassmanLineMinimization<T>::ptr2ions_);
    kbpsi_1.computeAll(*GrassmanLineMinimization<T>::ptr2ions_, orbitals1);

    KBPsiMatrixSparse kbpsi_2(
        GrassmanLineMinimization<T>::hamiltonian_->lapOper());
    kbpsi_2.setup(*GrassmanLineMinimization<T>::ptr2ions_);
    kbpsi_2.computeAll(*GrassmanLineMinimization<T>::ptr2ions_, orbitals2);

    // compute P^T*H*Q (orbitals1=P; orbitals2=Q)
    GrassmanLineMinimization<T>::mgmol_strategy_->computeHij(orbitals1,
        orbitals2, *GrassmanLineMinimization<T>::ptr2ions_, &kbpsi_1, &kbpsi_2,
        mat, consolidate);
}

// Compute P^T*H*P for orbitals1-->P and return result in mat.
// consolidate flag with either gather data or else return local partial
// contributions
template <class T>
void GrassmanCGSparse<T>::computeOrbitalsProdWithH(
    T& orbitals, VariableSizeMatrix<sparserow>& mat, const bool consolidate)
{
    // initialize KBPsiMatrices
    KBPsiMatrixSparse kbpsi(
        GrassmanLineMinimization<T>::hamiltonian_->lapOper());
    kbpsi.setup(*GrassmanLineMinimization<T>::ptr2ions_);
    kbpsi.computeAll(*GrassmanLineMinimization<T>::ptr2ions_, orbitals);

    // compute P^T*H*Q (orbitals1=P; orbitals2=Q)
    GrassmanLineMinimization<T>::mgmol_strategy_->computeHij(orbitals, orbitals,
        *GrassmanLineMinimization<T>::ptr2ions_, &kbpsi, mat, consolidate);

    return;
}

// parallel transport of history data
// update G=grad_, MG=pcgrad_ and Zo=sdir_
template <class T>
void GrassmanCGSparse<T>::parallelTransportUpdate(
    const double /*lambda*/, T& /*orbitals*/)
{
    // update gradient information
    //    SquareLocalMatrices<MATDTYPE,MemorySpace::Host> ss(grad->subdivx(),
    //    grad->chromatic_number()); grad->getLocalOverlap(ss);
    //    proj_matrices_->applyInvS(ss);
}

template class GrassmanCGSparse<LocGridOrbitals<ORBDTYPE>>;
