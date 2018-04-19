// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef KBPSIMATRIX_INTERFACE_H
#define KBPSIMATRIX_INTERFACE_H

#include "SparseDistMatrix.h"
#include "RemoteTasksDistMatrix.h"
#include "VariableSizeMatrix.h"
#include "tools.h"
#include "Timer.h"
#include "LocGridOrbitals.h"

//class LocGridOrbitals;
class ProjectedMatrices;
class ProjectedMatricesInterface;
class Ions;
class Ion;
class ProjectedMatricesSparseAOMM;

class KBPsiMatrixInterface
{
    int iterative_index_;
    
protected:
    static Timer   computeLocalElement_tm_;

public:
    KBPsiMatrixInterface()
        :iterative_index_(-1)
    {};
    
    virtual ~KBPsiMatrixInterface(){};
    
    int getIterativeIndex()const{ return iterative_index_;}
    void setOutdated(){ iterative_index_=-1; }
    void setIterativeIndex(const int index)
    {
        assert( index>=0 );
#ifdef DEBUG
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        iterative_index_=index;
    }
    
    virtual void addKBPsi(const int gid, const int st, const double val)=0;
    virtual void addKBBPsi(const int gid, const int st, const double val)=0;
    
    virtual double getEvnl(const Ions& ions, LocGridOrbitals& orbitals, ProjectedMatricesInterface* proj_matrices)=0;
    virtual double getValIonState(const int gid, const int st)const=0;
    virtual void scaleWithKBcoeff(const Ions& ions)=0;

    virtual void computeHvnlMatrix(const Ions&,dist_matrix::SparseDistMatrix<DISTMATDTYPE>&)const=0;
    virtual void computeHvnlMatrix(const Ions&,ProjectedMatricesInterface*)const=0;

    virtual void computeAll(Ions& ions, LocGridOrbitals& orbitals)=0;
    virtual void setup(const Ions& ions, 
               const LocGridOrbitals& orbitals)=0;

    /* Default implementation - returns error message */
    virtual void registerRemoteTasksDistMatrix(dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>* remote_tasks_DistMatrix)
    {
        (void) remote_tasks_DistMatrix;
        

        exitWithErrorMessage("registerRemoteTasksDistMatrix");
    }     
    virtual void computeHvnlMatrix(const Ions& ions,VariableSizeMatrix<sparserow>& vsmat)const
    {
        (void) ions;
        (void) vsmat;
        
        exitWithErrorMessage("registerRemoteTasksDistMatrix");        
    }        
    
    virtual void printTimers(ostream& os);

    void computeLocalElement(Ion& ion, const int istate, 
                     const int iloc, 
                     const ORBDTYPE* const psi, 
                     const bool flag);
};

#endif
