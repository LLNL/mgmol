// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "Hamiltonian.h"
#include "Control.h"
#include "Mesh.h"
#include "LocGridOrbitals.h"
#include "Potentials.h"
#include "ProjectedMatrices.h"

Timer Hamiltonian::apply_Hloc_tm_("Hamiltonian::apply_Hloc");

Hamiltonian::Hamiltonian()
{
    itindex_=-1;
    lapOper_ = NULL;
    hlphi_   = NULL;
    pot_     = new Potentials();
};

Hamiltonian::~Hamiltonian()
{
    if( hlphi_ != NULL )delete hlphi_;
    if( lapOper_ != NULL )delete lapOper_;
    delete pot_;
}

void Hamiltonian::setup(const pb::Grid& myGrid, const int lap_type)
{
    if( lapOper_ != NULL )delete lapOper_;
    lapOper_ = LapFactory<ORBDTYPE>::createLap(myGrid,lap_type);
}

const LocGridOrbitals& Hamiltonian::applyLocal(LocGridOrbitals& phi, const bool force)
{
    assert( phi.getIterativeIndex()>=0 );
    assert( pot_->getIterativeIndex()>=0 );

    if( hlphi_==NULL )
        hlphi_=new LocGridOrbitals(phi,false);
    if( !hlphi_->isCompatibleWith(phi) )
    {
        delete hlphi_;
        itindex_=-1;
        hlphi_=new LocGridOrbitals(phi,false);
    }
    const int new_index=100*phi.getIterativeIndex()+pot_->getIterativeIndex();
#ifdef DEBUG
    if( onpe0 ){
        (*MPIdata::sout)<<"Hamiltonian::applyLocal(), new_index ="<<new_index<<endl;
        (*MPIdata::sout)<<"Hamiltonian::applyLocal(), itindex_  ="<<itindex_<<endl;
    }
#endif
    if( force || new_index!=itindex_)
    {
        applyLocal(0, phi.chromatic_number(), phi, *hlphi_ );
    
        itindex_=new_index;
#ifdef PRINT_OPERATIONS
    }else{
        if( onpe0 )
            (*MPIdata::sout)<<"Hamiltonian::hlphi up to date, itindex_="<<itindex_<<endl;    
#endif
    }
    return *hlphi_;
}

void Hamiltonian::applyLocal(const int first_state, 
                             const int ncolors,
                             LocGridOrbitals& phi, 
                             LocGridOrbitals& hphi )
{
    apply_Hloc_tm_.start();
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        (*MPIdata::sout)<<"Hamiltonian::applyLocal() for states "<<first_state
            <<" to "<<first_state+ncolors-1<<endl;    
#endif
    assert(first_state>-1);
    
    const Control& ct = *(Control::instance());
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid  = mymesh->grid();

    const POTDTYPE* const vtot=pot_->vtot();

    phi.setDataWithGhosts();
    phi.trade_boundaries();
    
    if( ct.Mehrstellen() )
    {
        pb::GridFunc<POTDTYPE>  gfpot(mygrid,ct.bc[0],ct.bc[1],ct.bc[2]);
        gfpot.assign(vtot);
        if( ct.Mehrstellen() ) gfpot.trade_boundaries();
        const vector<vector<int> >& gid( phi.getGlobalIndexes() );
        pb::GridFuncVector<ORBDTYPE> gfvw1(false,mygrid,ct.bc[0],ct.bc[1],ct.bc[2],gid);
        pb::GridFuncVector<ORBDTYPE> gfvw2(false,mygrid,ct.bc[0],ct.bc[1],ct.bc[2],gid);
        //if( onpe0 )(*MPIdata::sout)<<"Hamiltonian::applyLocal, index="<<phi.getIterativeIndex()<<endl;
        for(int i=0;i<ncolors;i++)
        {
            pb::GridFunc<ORBDTYPE>*  gfw=new pb::GridFunc<ORBDTYPE>(mygrid,ct.bc[0],ct.bc[1],ct.bc[2]);
            gfvw1.push_back(gfw);
            gfvw2.push_back(&phi.getFuncWithGhosts(first_state+i));
        }
        gfvw1.prod(gfvw2,gfpot);
    
        pb::GridFunc<ORBDTYPE>  gf_work1(mygrid,ct.bc[0],ct.bc[1],ct.bc[2]);
        pb::GridFunc<ORBDTYPE>  gf_work2(mygrid,ct.bc[0],ct.bc[1],ct.bc[2]);
        for(int i=0;i<ncolors;i++)
        {
            // work1 = B*V*psi
            lapOper_->rhs(gfvw1.func(i),gf_work1);

            // work2 = -Lap*phi
            lapOper_->apply(phi.getFuncWithGhosts(first_state+i),
                            gf_work2);

            gf_work1+=gf_work2;
            hphi.setPsi(gf_work1,i+first_state);       
        }
        for(int i=0;i<ncolors;i++)
        {
            delete &gfvw1.func(i);
        }

    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for  
#endif
        for(int i=0;i<ncolors;i++)
        {
            // work1 = -Lap*phi
            lapOper_->applyWithPot(phi.getFuncWithGhosts(first_state+i),vtot,
                                   hphi.getPsi(i+first_state));
        }

    }
    
    apply_Hloc_tm_.stop();
}

// add to hij the elements <phi1|Hloc|phi2>
// corresponding to the local part of the Hamiltonian
void Hamiltonian::addHlocal2matrix(LocGridOrbitals& phi1, 
                              LocGridOrbitals& phi2, 
                              dist_matrix::SparseDistMatrix<DISTMATDTYPE>& hij, const bool force)
{
    applyLocal(phi2, force);

#ifdef PRINT_OPERATIONS
    if( onpe0 )
        (*MPIdata::sout)<<"Hamiltonian::addHlocal2matrix()"<<endl;
#endif
 
    phi1.addDotWithNcol2Matrix(*hlphi_, hij);
}

void Hamiltonian::addHlocalij(LocGridOrbitals& phi1, 
                              LocGridOrbitals& phi2)
{
    applyLocal(phi2);

#ifdef PRINT_OPERATIONS
    if( onpe0 )
        (*MPIdata::sout)<<"Hamiltonian::addHLocalij()"<<endl;
#endif
 
    phi1.addDot2H(*hlphi_);
}

void Hamiltonian::addHlocal2matrix(LocGridOrbitals& phi1, 
                              LocGridOrbitals& phi2,
                              VariableSizeMatrix<sparserow>& mat, const bool force)
{
    Control& ct = *(Control::instance());
 
    applyLocal(phi2, force);

#ifdef PRINT_OPERATIONS
    if( onpe0 )
        (*MPIdata::sout)<<"Hamiltonian::addHLocalij()"<<endl;
#endif

    SquareLocalMatrices<MATDTYPE> ss(phi1.subdivx(), phi1.chromatic_number());
 
    phi1.addDot2H(*hlphi_, ss);
    mat.initializeMatrixElements(ss, phi1.getGlobalIndexes(), ct.numst);
}

