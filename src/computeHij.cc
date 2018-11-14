// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cassert>

#include "MGmol_blas1.h"
#include "DistMatrix.h"

#include "Lap.h"
#include "LocGridOrbitals.h"
#include "ProjectedMatricesInterface.h"
#include "Ions.h"
#include "Hamiltonian.h"
#include "Potentials.h"
#include "ProjectedMatricesSparse.h"
#include "KBPsiMatrixSparse.h"
#include "GridFuncVector.h"
#include "Control.h"
#include "Mesh.h"
#include "KBPsiMatrix.h"
#include "Energy.h"
#include "MGmol.h"
#include "MGmol_MPI.h"
#include "DistMatrixWithSparseComponent.h"


static int sparse_distmatrix_nb_partitions=128;

void MGmol::computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrixSparse* const kbpsi_i,
            const KBPsiMatrixSparse* const kbpsi_j,
            VariableSizeMatrix<sparserow>& mat, const bool consolidate)
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        os_<<"computeHij()"<<endl;
#endif

    //compute phi_i^T*Hnl*Phi_j
    kbpsi_i->computeHvnlMatrix(kbpsi_j,ions,mat);

    // add local Hamiltonian part to phi_i^T*H*phi_j 
    addHlocal2matrix(orbitals_i,
                orbitals_j, 
                mat);

    // sum matrix elements among processors
    if(consolidate)
    {
        //get data distribution object
        DataDistribution *distributorH = BasicDataDistributors::orbitalsProdWithHDistributor();
        //gather data for only centered data
        distributorH->augmentLocalData(mat, false);
        
        //get indexes of centered data
        std::vector<int>locfcns;
        (*lrs_).getLocalSubdomainIndices(locfcns);
        
        //sparsify matrix rows according to pattern/ clear rows corresponding to non-centered data
        vector<int>pattern(mat.n(),0);
        for(std::vector<int>::iterator it=locfcns.begin(); it!=locfcns.end(); ++it)
        {
           const int *rindex = (int *)mat.getTableValue(*it);
           pattern[*rindex] = 1;
        }
        mat.sparsify(pattern);
        
        //gather/ distribute data from neighbors whose centered functions 
        //overlap with functions centered on local subdomain
        Mesh* mymesh = Mesh::instance();
        const pb::Grid& mygrid  = mymesh->grid();  
        const pb::PEenv& myPEenv=mymesh->peenv();
        double domain[3]={mygrid.ll(0),mygrid.ll(1),mygrid.ll(2)};
        DataDistribution distributor("Hij",2*(*lrs_).max_radii(), myPEenv, domain); 
        distributor.updateLocalRows(mat, true);  
    }
}

void MGmol::computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrixSparse* const kbpsi,
            VariableSizeMatrix<sparserow>& mat, const bool consolidate)
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        os_<<"computeHij()"<<endl;
#endif

    kbpsi->computeHvnlMatrix(ions,mat);

    // add local Hamiltonian part to phi^T*H*phi 
    addHlocal2matrix(orbitals_i,
                orbitals_j, 
                mat);

    // sum matrix elements among processors
    if(consolidate)
    {
        //get data distribution object
        DataDistribution *distributorH = BasicDataDistributors::orbitalsProdWithHDistributor();
        //gather data for only centered data
        distributorH->augmentLocalData(mat, false);
        
        //get indexes of centered data
        std::vector<int>locfcns;
        (*lrs_).getLocalSubdomainIndices(locfcns);
        
        //sparsify matrix rows according to pattern/ clear rows corresponding to non-centered data
        vector<int>pattern(mat.n(),0);
        for(std::vector<int>::iterator it=locfcns.begin(); it!=locfcns.end(); ++it)
        {
           const int *rindex = (int *)mat.getTableValue(*it);
           pattern[*rindex] = 1;
        }
        mat.sparsify(pattern);
        
        //gather/ distribute data from neighbors whose centered functions 
        //overlap with functions centered on local subdomain
        Mesh* mymesh = Mesh::instance();
        const pb::Grid& mygrid  = mymesh->grid();  
        const pb::PEenv& myPEenv=mymesh->peenv();
        double domain[3]={mygrid.ll(0),mygrid.ll(1),mygrid.ll(2)};
        DataDistribution distributor("Hij",2*(*lrs_).max_radii(), myPEenv, domain); 
        distributor.updateLocalRows(mat, true);  
    }
}

void MGmol::computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrixSparse* const kbpsi_i,
            const KBPsiMatrixSparse* const kbpsi_j,
            dist_matrix::DistMatrix<DISTMATDTYPE>& hij)
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        os_<<"computeHij()"<<endl;
#endif
#ifdef USE_MPI
    assert( remote_tasks_DistMatrix_!=0 );
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm=mmpi.commSameSpin();
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sparseH(comm,hij,
                                     remote_tasks_DistMatrix_,
                                     sparse_distmatrix_nb_partitions);
#else
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sparseH(0,hij,
                                     rtasks_distmatrix);
#endif

    kbpsi_i->computeHvnlMatrix(kbpsi_j,ions,sparseH);

    // add local Hamiltonian part to phi^T*H*phi 
    addHlocal2matrix(orbitals_i,
                orbitals_j, 
                sparseH);

    // sum matrix elements among processors
    sparseH.parallelSumToDistMatrix();
}

void MGmol::computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrix* const kbpsi_i,
            const KBPsiMatrix* const kbpsi_j,
            dist_matrix::DistMatrix<DISTMATDTYPE>& hij)
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        os_<<"computeHij()"<<endl;
#endif
#ifdef USE_MPI
    assert( remote_tasks_DistMatrix_!=0 );
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm=mmpi.commSameSpin();
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sparseH(comm,hij,
                                     remote_tasks_DistMatrix_,
                                     sparse_distmatrix_nb_partitions);
#else
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sparseH(0,hij,
                                     rtasks_distmatrix);
#endif

    kbpsi_i->computeHvnlMatrix(kbpsi_j,ions,sparseH);

    // add local Hamiltonian part to phi^T*H*phi 
    addHlocal2matrix(orbitals_i,
                orbitals_j, 
                sparseH);

    // sum matrix elements among processors
    sparseH.parallelSumToDistMatrix();
}

void MGmol::computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrixInterface* const kbpsi,
            dist_matrix::DistMatrix<DISTMATDTYPE>& hij)
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        os_<<"computeHij()"<<endl;
#endif
#ifdef USE_MPI
    assert( remote_tasks_DistMatrix_!=0 );
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm=mmpi.commSameSpin();
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sparseH(comm,hij,
                                     remote_tasks_DistMatrix_,
                                     sparse_distmatrix_nb_partitions);
#else
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sparseH(0,hij,
                                     rtasks_distmatrix);
#endif

    kbpsi->computeHvnlMatrix(ions,sparseH);

    // add local Hamiltonian part to phi^T*H*phi 
    addHlocal2matrix(orbitals_i,
                orbitals_j, 
                sparseH);

    // sum matrix elements among processors
    sparseH.parallelSumToDistMatrix();
}

void MGmol::computeHij(LocGridOrbitals& orbitals_i,
                      LocGridOrbitals& orbitals_j,
                      const Ions& ions,
                      const KBPsiMatrixInterface* const kbpsi,
                      dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sparseH)
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        os_<<"computeHij()"<<endl;
#endif

    kbpsi->computeHvnlMatrix(ions,sparseH);

    // add local Hamiltonian part to phi^T*H*phi 
    addHlocal2matrix(orbitals_i,
                orbitals_j, 
                sparseH);

    // sum matrix elements among processors
    sparseH.parallelSumToDistMatrix();
}

void MGmol::computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrixInterface* const kbpsi,
            ProjectedMatricesInterface* projmatrices)
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        os_<<"computeHij()"<<endl;
#endif

    kbpsi->computeHvnlMatrix(ions,projmatrices);

    // add local Hamiltonian part to phi^T*H*phi 
    addHlocalij(orbitals_i,
                orbitals_j, projmatrices);

    projmatrices->consolidateH();
}

void MGmol::getKBPsiAndHij(LocGridOrbitals& orbitals_i,
                    LocGridOrbitals& orbitals_j,
                    Ions& ions,
                    KBPsiMatrix* kbpsi_i,
                    KBPsiMatrix* kbpsi_j,
                    ProjectedMatricesInterface* projmatrices,
                    dist_matrix::DistMatrix<DISTMATDTYPE>& hij)
{    
    kbpsi_i->computeAll(ions, orbitals_i);
    if( kbpsi_j != kbpsi_i )
        kbpsi_j->computeAll(ions, orbitals_j);

    computeHij(orbitals_i, orbitals_j, ions, 
           kbpsi_i, kbpsi_j, hij);

    projmatrices->setHiterativeIndex(orbitals_j.getIterativeIndex(), 
                                    hamiltonian_->potential().getIterativeIndex());
}

void MGmol::getKBPsiAndHij(LocGridOrbitals& orbitals_i,
                    LocGridOrbitals& orbitals_j,
                    Ions& ions,
                    KBPsiMatrixInterface* kbpsi,
                    ProjectedMatricesInterface* projmatrices,
                    dist_matrix::DistMatrix<DISTMATDTYPE>& hij)
{    
    kbpsi->computeAll(ions, orbitals_i);

    computeHij(orbitals_i, orbitals_j, ions, kbpsi, hij);

    projmatrices->setHiterativeIndex(orbitals_j.getIterativeIndex(), 
                                    hamiltonian_->potential().getIterativeIndex());
}

void MGmol::getKBPsiAndHij(LocGridOrbitals& orbitals_i,
                    LocGridOrbitals& orbitals_j,
                    Ions& ions,
                    KBPsiMatrixInterface* kbpsi,
                    ProjectedMatricesInterface* projmatrices,
                    dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sh)
{    
    kbpsi->computeAll(ions, orbitals_i);

    sh.clearData();
    computeHij(orbitals_i, orbitals_j, ions, kbpsi, sh);

    projmatrices->setHiterativeIndex(orbitals_j.getIterativeIndex(), 
                                    hamiltonian_->potential().getIterativeIndex());
}

void MGmol::getKBPsiAndHij(LocGridOrbitals& orbitals_i,
                    LocGridOrbitals& orbitals_j,
                    Ions& ions,
                    KBPsiMatrixInterface* kbpsi,
                    ProjectedMatricesInterface* projmatrices)
{    
    kbpsi->computeAll(ions, orbitals_i);

    projmatrices->clearSparseH();

    computeHij(orbitals_i, orbitals_j, ions, kbpsi, projmatrices);

    projmatrices->setHiterativeIndex(orbitals_j.getIterativeIndex(), 
                                    hamiltonian_->potential().getIterativeIndex());
}

void MGmol::getKBPsiAndHij(LocGridOrbitals& orbitals,
                    Ions& ions,
                    KBPsiMatrixInterface* kbpsi,
                    dist_matrix::DistMatrix<DISTMATDTYPE>& hij)
{
    getKBPsiAndHij(orbitals, orbitals, ions,
                   kbpsi,
                   proj_matrices_, hij);
}

void MGmol::getKBPsiAndHij(LocGridOrbitals& orbitals,
                           Ions& ions)
{
    getKBPsiAndHij(orbitals, orbitals, ions,
                   g_kbpsi_,
                   proj_matrices_);
}

void MGmol::computeHnlPhiAndAdd2HPhi(Ions& ions,
                LocGridOrbitals& phi,
                LocGridOrbitals& hphi,
                const KBPsiMatrixInterface* const kbpsi) 
{
    // H_nl
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        os_<<"computeHnlPhiAndAdd2HPhi()"<<endl;
#endif

    Control& ct = *(Control::instance());
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();
    const int numpt=mygrid.size();
    
    const vector<vector<int> >& gid( phi.getOverlappingGids() );
    const short ncolors=gid[0].size();
 
    {
        // compute Hnl*phi
        if( ct.Mehrstellen() )
        {
            ORBDTYPE* work1=new ORBDTYPE[numpt];
            pb::GridFuncVector<ORBDTYPE> gfv(false,mygrid,ct.bc[0],ct.bc[1],ct.bc[2],gid);
            for(short icolor=0;icolor<ncolors;icolor++)
            {
                get_vnlpsi(ions, gid, icolor, kbpsi, work1);
                pb::GridFunc<ORBDTYPE>*  gf_work
                    =new pb::GridFunc<ORBDTYPE>(mygrid,ct.bc[0],ct.bc[1],ct.bc[2]);
                gfv.push_back(gf_work);
                gf_work->assign(work1);
            }
            delete[] work1;
            gfv.trade_boundaries();
            
            // compute B*Hnl*phi and add it to H*phi
            ORBDTYPE* hnl=new ORBDTYPE[numpt];
            for(int icolor=0;icolor<ncolors;icolor++)
            {
                hamiltonian_->lapOper()->rhs(gfv.func(icolor),hnl);
 
                // Add the contribution of the non-local potential to H phi
                ORBDTYPE* hpsi=hphi.getPsi(icolor);
                MPaxpy(numpt, 1., hnl, hpsi);
            }
            delete[] hnl;
 
            for(short color=0;color<ncolors;color++)
                delete &gfv.func(color);
        }
        else // no Mehrstellen
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
{
            ORBDTYPE* hnl=new ORBDTYPE[numpt];
#ifdef _OPENMP
#pragma omp for
#endif
            for(short icolor=0;icolor<ncolors;icolor++)
            {
                get_vnlpsi(ions, gid, icolor, kbpsi, hnl);
                MPaxpy(numpt, 1., hnl, hphi.getPsi(icolor));
            }
            
            delete[] hnl;
}
        }

    }

    hphi.setIterativeIndex( phi.getIterativeIndex() );
}

void MGmol::addHlocal2matrix(LocGridOrbitals& orbitalsi,
                             LocGridOrbitals& orbitalsj, 
                             VariableSizeMatrix<sparserow>& mat)
{
    computeHij_tm_.start();
    
#if DEBUG
    os_<<" addHlocal2matrix()"<<endl;
#endif

    hamiltonian_->addHlocal2matrix(orbitalsi, orbitalsj, mat, true);
    
    computeHij_tm_.stop();
}

void MGmol::addHlocal2matrix(
                LocGridOrbitals& orbitalsi,
                LocGridOrbitals& orbitalsj, 
                dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sparseH)
{
    computeHij_tm_.start();
    
#if DEBUG
    os_<<" addHlocal2matrix()"<<endl;
#endif

    hamiltonian_->addHlocal2matrix(orbitalsi, orbitalsj, sparseH);
    
    computeHij_tm_.stop();
}

void MGmol::addHlocal2matrix(LocGridOrbitals& orbitalsi,
                        LocGridOrbitals& orbitalsj, 
                        dist_matrix::DistMatrixWithSparseComponent<DISTMATDTYPE>& mat)
{
    computeHij_tm_.start();
    
#if DEBUG
    os_<<" addHlocal2matrix()"<<endl;
#endif

    //create a new sparse data structure initialized with content of mat.sparse()
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sparse(mat.sparse());

    //add local H to sparse
    hamiltonian_->addHlocal2matrix(orbitalsi, orbitalsj, sparse);
    
    //sum up values in sparse into a DistMatrix (initialized with zeroes)
    sparse.parallelSumToDistMatrix();
    
    computeHij_tm_.stop();
}

void MGmol::addHlocalij(LocGridOrbitals& orbitalsi,
                        LocGridOrbitals& orbitalsj,
                        ProjectedMatricesInterface* pmat)
{
    computeHij_tm_.start();
    
#if DEBUG
    os_<<" addHlocalij()"<<endl;
#endif

    hamiltonian_->addHlocalij(orbitalsi, orbitalsj, pmat);
    
    computeHij_tm_.stop();
}

void MGmol::getHpsiAndTheta(Ions& ions, 
                      LocGridOrbitals& phi, 
                      LocGridOrbitals& hphi)
{
    getHpsiAndTheta(ions,phi,hphi,g_kbpsi_);
}

void MGmol::getHpsiAndTheta(Ions& ions, 
                      LocGridOrbitals& phi, 
                      LocGridOrbitals& hphi,
                      const KBPsiMatrixInterface* const kbpsi)
{
    get_Hpsi_and_Hij_tm_.start();
    
    const Potentials& pot =hamiltonian_->potential();
    
    const int phi_it_index=phi.getIterativeIndex();

#if DEBUG
    os_<<" getHpsiAndTheta"<<endl;
#endif

    hphi.assign( hamiltonian_->applyLocal(phi) );
 
    // Compute "nstates" columns of matrix
    //  Hij = phi**T * H_loc * phi  and save in sh
    if( proj_matrices_->isHupToDate(phi.getIterativeIndex(),pot.getIterativeIndex()) )
    {
#ifdef PRINT_OPERATIONS
        if( onpe0 )
            os_<<"Hij matrix up to date, no computation necessary"<<endl;
#endif
    }else
    {
#ifdef PRINT_OPERATIONS
        if( onpe0 )
            os_<<"build matrix Hij = Phi**T * H * Phi"<<endl;
#endif

        proj_matrices_->clearSparseH();

        // Compute the contribution of the non-local potential into sh
        kbpsi->computeHvnlMatrix(ions,proj_matrices_);

        // add local part of H to sh
        SquareLocalMatrices<MATDTYPE> slh(phi.subdivx(),phi.chromatic_number());

        phi.computeLocalProduct(hphi, slh);

        proj_matrices_->addMatrixElementsSparseH(slh); 

        proj_matrices_->consolidateH();
        
        energy_->saveVofRho();

        // Compute matrix theta = invB * Hij
        proj_matrices_->updateTheta();

        proj_matrices_->setHiterativeIndex(phi_it_index, pot.getIterativeIndex());
    }
    computeHnlPhiAndAdd2HPhi(ions,phi,hphi,kbpsi);
    hphi.setIterativeIndex(phi_it_index);

    get_Hpsi_and_Hij_tm_.stop();
}
