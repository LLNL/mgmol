// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "EigenDMStrategy.h"
#include "LocGridOrbitals.h"
#include "Control.h"
#include "ProjectedMatrices.h"

EigenDMStrategy::EigenDMStrategy(
    LocGridOrbitals* current_orbitals,
    ProjectedMatricesInterface* proj_matrices):
        current_orbitals_(current_orbitals),
        proj_matrices_(proj_matrices)
{
}

void EigenDMStrategy::initialize()
{
    update(); 
}

int EigenDMStrategy::update()
{
    Control& ct = *(Control::instance());

    dist_matrix::DistMatrix<DISTMATDTYPE>  zz("Z", ct.numst, ct.numst);
    
    ProjectedMatrices* pmat=dynamic_cast<ProjectedMatrices*>( proj_matrices_ );
    pmat->updateDMwithEigenstatesAndRotate(current_orbitals_->getIterativeIndex(),zz); 
    
    //if( onpe0 && ct.verbose>2 )
    //    (*MPIdata::sout)<<"get_dm_diag: rotate orbitals "<<endl;
    current_orbitals_->multiply_by_matrix(zz);
    current_orbitals_->setDataWithGhosts();
    current_orbitals_->trade_boundaries();
    
    return 0;
}
