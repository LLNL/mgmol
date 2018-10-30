// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifdef PROCRUSTES
#include "munkres.h"


void procrustes(dist_matrix::DistMatrix<DISTMATDTYPE>& a,
                dist_matrix::DistMatrix<DISTMATDTYPE>& b,
                dist_matrix::DistMatrix<DISTMATDTYPE>& p)
{
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = mbc. bcxt();
    const int nst    =a.m();

    p.gemm('t','n',1.,b,a,0.);
    
    munkres::Matrix<double> matrix(nst, nst);
    
    double maxval=0.;
    for ( int row = 0 ; row < nst ; row++ ) {
        for ( int col = 0 ; col < nst ; col++ ) {
            maxval = max( maxval, p.val(row+nst*col) );
        }
    }
    
    // Initialize matrix
    for ( int row = 0 ; row < nst ; row++ ) {
        for ( int col = 0 ; col < nst ; col++ ) {
            matrix(row,col) = maxval-fabs( p.val(row+nst*col) );
        }
    }

    // Apply Munkres algorithm to matrix.
    munkres::Munkres m;
    m.solve(matrix);

    // 0 -> 1, -1 -> 0
    for ( int row = 0 ; row < nst ; row++ ) {
        for ( int col = 0 ; col < nst ; col++ ) {
            matrix(row,col)+=1.;
        }
    }
    
    for ( int row = 0 ; row < nst ; row++ ) {
        for ( int col = 0 ; col < nst ; col++ ) {
            if( matrix(row,col)>0.5 ){
                if( p.val(row+nst*col)<0. )
                    matrix(row,col)=-1.;
                break;
            }
        }
    }
    
    for ( int row = 0 ; row < nst ; row++ ) {
        for ( int col = 0 ; col < nst ; col++ ) {
            p.setval(row+nst*col, matrix(row,col));
        }
    }

    //(*MPIdata::sout)<<"Procrustes matrix..."<<endl;
    //p.print((*MPIdata::sout),0,0,nst,nst); 

    //w.gemm('n','n',1.,b,p,0.);
}
#endif

