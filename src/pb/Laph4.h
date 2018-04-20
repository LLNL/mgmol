// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: Laph4.h,v 1.10 2010/01/28 22:56:47 jeanluc Exp $
#ifndef LAPH4_H
#define LAPH4_H

#include "Laph2.h"

// Laplacian operator O(h^4)
namespace pb{
template <class T>
class Laph4:public Lap<T>{
    double  diagEl_;
    double  invDiagEl_;

    Laph2<T>*  lower_order_op_;

public:

    Laph4(const Grid& mygrid):Lap<T>(mygrid){
        //cout<<" Create Laph4 operator\n";

        diagEl_   = 2.5*(Lap<T>::inv_h2(0)+Lap<T>::inv_h2(1)+Lap<T>::inv_h2(2)); // 2.5 = 30./12.
        invDiagEl_= 1./diagEl_;
        Lap<T>::name_="Classical 4th order";

        lower_order_op_=NULL;
    }
    
    ~Laph4()
    {
        if( lower_order_op_!=NULL ){
            delete lower_order_op_;
            lower_order_op_=NULL;
        }
    }

    // construct a coarse grid operator
    Laph4 coarseOp(const Grid& mygrid)
    {
        Grid coarse_G=mygrid.coarse_grid();

        Laph4 A(coarse_G);

        return A;
    }

    Laph4 replicatedOp(const Grid& replicated_grid)
    {
        Laph4 replicated_A(replicated_grid);
        
        return replicated_A;
    }

    void setLowerOrderGrid()
    {
        this->setFDLowerOrderGrid(Laph2<T>::minNumberGhosts());
    }
    
    Laph2<T>& getLowerOrderOp()
    {
        if( lower_order_op_==NULL ){
            this->setFDLowerOrderGrid(Laph2<T>::minNumberGhosts());
            lower_order_op_=new Laph2<T>( Lap<T>::getLowerOrderGrid() );
        }
        return *lower_order_op_;
    }

    static short minNumberGhosts()
    {
        return 2;
    }

    // A->B
    void apply(GridFunc<T> &A, GridFunc<T> &B)
    {
        this->del2_4th(A,B);
        B.set_bc(A.bc(0),A.bc(1),A.bc(2));
    }
    void applyWithPot(GridFunc<T> &A, const double* pot, T* B)
    {
        this->del2_4th_withPot(A,pot,B);
    }
    void apply(GridFuncVector<T>& A, GridFuncVector<T> &B)
    {
        assert( A.size()==B.size() );
        A.trade_boundaries();
        const int nfunc=(int)A.size();
        for(int k=0;k<nfunc;k++){
            this->del2_4th(A.func(k),B.func(k));
        }
    }
 
    void jacobi(GridFunc<T>&, const GridFunc<T>&, GridFunc<T>&);
    void jacobi(GridFuncVector<T>&, const GridFuncVector<T>&, GridFunc<T>&);
    void jacobi(GridFuncVector<T>&, const GridFuncVector<T>&, GridFuncVector<T>&);

    double diagEl(void)const{ return diagEl_; }; 
    double invDiagEl(void)const{ return invDiagEl_; }; 
};

} // namespace pb

#endif
