// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Control.h"
#include "HamiltonianMVP_DMStrategy.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesSparse.h"
#include "LocGridOrbitals.h"
#include "HamiltonianMVPSolver.h"
#include "Ions.h"
#include "DistMatrix.h"
#include "DistMatrixWithSparseComponent.h"
#include <vector>
using namespace std;

template <class T1, class T2, class T3>
HamiltonianMVP_DMStrategy<T1, T2, T3>::HamiltonianMVP_DMStrategy(
    MPI_Comm comm, 
    ostream& os, 
    Ions& ions,
    Rho* rho,
    Energy* energy,
    Electrostatic* electrostat,
    MGmol* mgmol_strategy,
    LocGridOrbitals* orbitals):
        comm_(comm),
        os_(os),
        ions_(ions),
        rho_(rho),
        energy_(energy),
        electrostat_(electrostat),
        global_indexes_(orbitals->getGlobalIndexes()),
        orbitals_(orbitals),
        mgmol_strategy_(mgmol_strategy)
{
    Control& ct = *(Control::instance());
    
    assert( electrostat_!=0 );
    assert( energy_!=0 );
    assert( ct.dm_inner_steps>0 );

    T3* projmatrices=dynamic_cast<T3*>( orbitals->getProjMatrices() );
    
    solver_ = new HamiltonianMVPSolver<T1, T2, T3>(comm_, os_,ions_,
             rho_, energy_, electrostat_,
             mgmol_strategy_,
             ct.numst,ct.occ_width,
             ct.getNel(),
             global_indexes_,
             ct.dm_inner_steps,
             projmatrices->getH(),
             true);
    
}

template <class T1, class T2, class T3>
HamiltonianMVP_DMStrategy<T1, T2, T3>::~HamiltonianMVP_DMStrategy()
{
    delete solver_;
}

template <class T1, class T2, class T3>
void HamiltonianMVP_DMStrategy<T1, T2, T3>::initialize()
{
}

template <class T1, class T2, class T3>
int HamiltonianMVP_DMStrategy<T1, T2, T3>::update()
{
    assert( solver_!=0 );
    
    Control& ct = *(Control::instance());
    if( onpe0 && ct.verbose>2 )
    {
        (*MPIdata::sout)<<"HamiltonianMVP_DMStrategy::update()..."<<endl;
    }

    return solver_->solve(*orbitals_);
}

template <class T1, class T2, class T3>
void HamiltonianMVP_DMStrategy<T1, T2, T3>::stripDM()
{
}

template <class T1, class T2, class T3>
void HamiltonianMVP_DMStrategy<T1, T2, T3>::dressDM()
{
}

template <class T1, class T2, class T3>
void HamiltonianMVP_DMStrategy<T1, T2, T3>::reset()
{
    solver_->reset();
}

template class HamiltonianMVP_DMStrategy< dist_matrix::DistMatrix<DISTMATDTYPE>, dist_matrix::DistMatrixWithSparseComponent<DISTMATDTYPE>, ProjectedMatrices >;
template class HamiltonianMVP_DMStrategy< VariableSizeMatrix<sparserow>, VariableSizeMatrix<sparserow>, ProjectedMatricesSparse >;