// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <iostream>
#include <iomanip>
using namespace std;

#include <string.h>
#include "DensityMatrix.h"
#include "MGmol_MPI.h"

const double factor_kernel4dot=10.;

#define PROCRUSTES 0

// occupations in [0,1]
// DM eigenvalues in [0,orbital_occupation]

DensityMatrix::DensityMatrix(const int ndim)
{
    assert( ndim>=0 );
    
    dim_=ndim;

    occ_uptodate_=false;
    stripped_    =false;
    uniform_occ_ = false;

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    orbital_occupation_= mmpi.nspin() >1 ? 1. : 2.;

    orbitals_index_=-1;

    if( dim_>0 )
    {
        dm_         =new dist_matrix::DistMatrix<DISTMATDTYPE>("DM",    ndim, ndim);
        kernel4dot_ =new dist_matrix::DistMatrix<DISTMATDTYPE>("K4dot", ndim, ndim);
        work_       =new dist_matrix::DistMatrix<DISTMATDTYPE>("work",  ndim, ndim);
        occupation_.resize(dim_);
    }else{
        occ_uptodate_=true;
    }
    
    //if( onpe0 )
    //    (*MPIdata::sout)<<"DensityMatrix: factor_kernel4dot="<<factor_kernel4dot<<endl;
}

DensityMatrix::DensityMatrix(const DensityMatrix& dm)
{
    dim_=dm.dim_;

    uniform_occ_ = dm.uniform_occ_;
    occ_uptodate_=dm.occ_uptodate_;
    stripped_    =dm.stripped_;

    orbitals_index_=dm.orbitals_index_;

    if( dim_>0 ){
        dm_         =new dist_matrix::DistMatrix<DISTMATDTYPE>(*dm.dm_);
        kernel4dot_ =new dist_matrix::DistMatrix<DISTMATDTYPE>(*dm.kernel4dot_);
        work_       =new dist_matrix::DistMatrix<DISTMATDTYPE>(*dm.work_);
        occupation_.resize(dim_);
    }
    
    //if( onpe0 )
    //    (*MPIdata::sout)<<"DensityMatrix: factor_kernel4dot="<<factor_kernel4dot<<endl;
}

DensityMatrix& DensityMatrix::operator=(const DensityMatrix& dm)
{
    if( this==&dm )return *this;
    
    dim_=dm.dim_;

    occ_uptodate_=dm.occ_uptodate_;
    stripped_    =dm.stripped_;
    uniform_occ_ = dm.uniform_occ_;

    orbitals_index_=dm.orbitals_index_;

    *dm_         = *dm.dm_;
    *kernel4dot_ = *dm.kernel4dot_;
    *work_       = *dm.work_;
    occupation_=dm.occupation_;
    
    return *this;
}

DensityMatrix::~DensityMatrix()
{
    if( dim_>0 ){
        assert( dm_!=0 );
        assert( kernel4dot_!=0 );
        assert( work_!=0 );
        
        delete dm_;
        delete kernel4dot_;
        delete work_;
    }
}

void DensityMatrix::build(const dist_matrix::DistMatrix<DISTMATDTYPE>& zmat,
                          const vector<DISTMATDTYPE>& occ,
                          const int new_orbitals_index)
{    
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        (*MPIdata::sout)<<"DensityMatrix::build(const DistMatrix<DISTMATDTYPE>&,const vector<DISTMATDTYPE>&,const int)"<<endl;
#endif
    
    //diagonal matrix with occ values in diagonal
    dist_matrix::DistMatrix<DISTMATDTYPE> gamma("Gamma", &occ[0], dim_, dim_);
    gamma.scal(orbital_occupation_); // rescale for spin

    // work_ = zmat*gamma with gamma symmetric
    work_->symm('r', 'l', 1., gamma, zmat, 0.);
    
    // dm_ = work_ * zmat^T
    dm_->gemm('n', 't', 1., *work_, zmat, 0.);
    
    vector<DISTMATDTYPE> w(dim_);
    for(int i=0;i<dim_;i++)
        w[i]=(DISTMATDTYPE)(orbital_occupation_*min(1.,factor_kernel4dot*occ[i]));
    gamma.setDiagonal(w);
    
    work_->symm('r', 'l', 1., gamma, zmat, 0.);
    kernel4dot_->gemm('n', 't', 1., *work_, zmat, 0.);
    
    stripped_    =false;
    orbitals_index_=new_orbitals_index;
}

void DensityMatrix::build(const dist_matrix::DistMatrix<DISTMATDTYPE>& zmat,
                          const int new_orbitals_index)
{
    build(zmat, occupation_, new_orbitals_index);
}

// build diagonal matrix
void DensityMatrix::build(const vector<DISTMATDTYPE>& occ,
                          const int new_orbitals_index)
{
    assert( dm_!=0 );
    
    setOccupations(occ);

    build(new_orbitals_index);
} 

void DensityMatrix::build(const int new_orbitals_index)
{
//#ifdef PRINT_OPERATIONS
    if( onpe0 )
        (*MPIdata::sout)<<"DensityMatrix::build() for diagonal occupation..."<<endl;
//#endif
    if( !occ_uptodate_ && onpe0 )
        (*MPIdata::sout)<<"Warning: occupations not up to date to build DM!!!"<<endl;

    dist_matrix::DistMatrix<DISTMATDTYPE> gamma("Gamma", &occupation_[0], dim_, dim_);
    gamma.scal(orbital_occupation_); // rescale for spin

    *dm_=gamma;
    
    kernel4dot_->clear();
    vector<DISTMATDTYPE> w(dim_);
    for(int i=0;i<dim_;i++)
        w[i]=(DISTMATDTYPE)(orbital_occupation_*min(1.,factor_kernel4dot*occupation_[i]));
    kernel4dot_->setDiagonal(w);
    
    stripped_    =false;
    orbitals_index_=new_orbitals_index;
}
void DensityMatrix::setUniform(const DISTMATDTYPE nel, const int new_orbitals_index)
{    
    const DISTMATDTYPE occ=(DISTMATDTYPE)((double)nel/(double)dim_);
    assert( occ<1.01 );
    for(int i=0;i<dim_;i++)
        occupation_[i]=occ;

    occ_uptodate_=true;
    
    uniform_occ_=true;

    build(occupation_,new_orbitals_index);    
}

void DensityMatrix::buildFromBlock(const dist_matrix::DistMatrix<DISTMATDTYPE>& block00)
{
    dm_->clear();
    dm_->assign(block00,0,0);
    dm_->print((*MPIdata::sout),0,0,25,25);
}

void DensityMatrix::rotate(const dist_matrix::DistMatrix<DISTMATDTYPE>&  rotation_matrix,
                           const bool flag_eigen)
{

    if( !flag_eigen ){
        dist_matrix::DistMatrix<DISTMATDTYPE> invU(rotation_matrix);
        vector<int> ipiv;
        invU.getrf(ipiv);
        
        // dm -> u**-1 * dm
        invU.getrs('n',*dm_,ipiv);
        
        // tmp = dm**T * u**-T
        dist_matrix::DistMatrix<DISTMATDTYPE> tmp(rotation_matrix);
        tmp.transpose(*dm_);
        
        // tmp = u**-1 * dm * u**-T
        invU.getrs('n',tmp,ipiv);
        
        *dm_=tmp;
    }

}

void DensityMatrix::printOccupations(ostream& os)const
{
    if( onpe0 ){
        os<<endl<<" Occupation numbers: ";

        // Print ten to a row.
        os.setf(ios::right,ios::adjustfield);
        os.setf(ios::fixed,ios::floatfield);
        os<<setprecision(3);
        for(int i=0;i<dim_;i++){
            if ( (i%10) == 0 ) os << endl;
            os<<setw(7)<<occupation_[i]*orbital_occupation_<<" ";
        }

        os<<endl;
    }
}

//double DensityMatrix::getSumOccupations()const
//{
//    double sum=0.;
//    for(int i=0;i<dim_;i++)
//    {
//        sum+=(double)occupation_[i]*(double)orbital_occupation_;
//    }
//
//    return sum;
//}

// solve the eigenvalue problem L^T*dm_*L*V=occ*V
// using the LL^T decomposition of S to get occ
void DensityMatrix::diagonalize(const dist_matrix::DistMatrix<DISTMATDTYPE>& ls,
                                vector<DISTMATDTYPE>& occ)
{
    dist_matrix::DistMatrix<DISTMATDTYPE> evect("EigVect", dim_, dim_);
    
    *work_ = (*dm_);

    // *work_ = L**T * *work_
    work_->trmm('l','l','t','n',1.,ls); 
    // *work_ = *work_*L
    work_->trmm('r','l','n','n',1.,ls);

    // compute eigenvalues of work_
    work_->syev('n', 'l', occ, evect);
}

void DensityMatrix::diagonalize(const char eigv,
                                vector<DISTMATDTYPE>& occ,
                                dist_matrix::DistMatrix<DISTMATDTYPE>& vect)
{
    *work_ = (*dm_);

    // compute eigenvectors and eigenvalues of work_
    work_->syev(eigv, 'l', occ, vect);
}


void DensityMatrix::computeOccupations(const dist_matrix::DistMatrix<DISTMATDTYPE>& ls)
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        (*MPIdata::sout)<<"DensityMatrix::computeOccupations()"<<endl;
#endif
    
    vector<DISTMATDTYPE> occ(dim_);
    
    diagonalize(ls, occ);
    const double occinv=1./orbital_occupation_;
    
    const double tol=1.e-5;
    bool flag=false;
    for(int i=0;i<dim_;i++)
    {
        double occ_val = (double)occ[i]*occinv;
        (*MPIdata::sout)<<setprecision(16);
        if( onpe0 
         && (occ_val<0.-tol || occ_val>1.+tol)
            ){
            (*MPIdata::sout)<<"WARNING: DensityMatrix::computeOccupations(), occ["<<i<<"]="<<occ_val;
            //if( occ_uptodate_)(*MPIdata::sout)<<" vs. "<<occupation_[dim_-i-1];
            (*MPIdata::sout)<<endl;
            flag=true;
        }
        assert( occ_val>0.-tol );
        assert( occ_val<1.+tol );
        occ[i]=(DISTMATDTYPE)max(0.,occ_val);
        occ[i]=(DISTMATDTYPE)min(1.,occ_val);
    }
    if( flag )printOccupations((*MPIdata::sout));

    for(int i=0;i<dim_;i++)
    {
        occupation_[i]=occ[dim_-i-1];
    } 
    occ_uptodate_=true;
}

void DensityMatrix::setOccupations(const vector<DISTMATDTYPE>& occ)
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        (*MPIdata::sout)<<"DensityMatrix::setOccupations()"<<endl;
#endif
    assert( (int)occ.size()==dim_ );
    memcpy(&occupation_[0],&occ[0],dim_*sizeof(DISTMATDTYPE));
    occ_uptodate_=true;
}
    
double DensityMatrix::computeEntropy()const
{
    double s=0.;
    const double tol=1.e-15;
    const double tol_interval=1.e-6;
 
    assert( occ_uptodate_ );
    for (int st = 0; st < dim_; st++)
    {
        const double fi = (double)occupation_[st];
        if( fi>1.+tol_interval )
            (*MPIdata::sout)<<setprecision(15)<<scientific<<"f["<<st<<"]="<<fi<<endl;
        assert( fi>=0.-tol_interval );
        assert( fi<=1.+tol_interval );
        if(fi<tol){
            s += (1.-fi)*log(1.-fi);
        }else if( fi>1.-tol ){
            s += fi*log(fi);
        }else{
            s += fi*log(fi)+(1.-fi)*log(1.-fi);
        }
    }
    
    return (double)(-orbital_occupation_)*s; // in units of kbt
}

void DensityMatrix::setto2InvS(const dist_matrix::DistMatrix<DISTMATDTYPE>& invS,
                               const int orbitals_index)
{
    *dm_=invS;
    dm_->scal(orbital_occupation_);

    if( !occ_uptodate_ )
    {
        for(int st = 0; st < dim_; st++)
            occupation_[st]=1.;
        occ_uptodate_=true;
    }
    uniform_occ_=false;
    orbitals_index_=orbitals_index; 
}

void DensityMatrix::stripS(const dist_matrix::DistMatrix<DISTMATDTYPE>& ls,
                           const int orbitals_index_gram)
{
    assert( !stripped_ );
    //assert( orbitals_index_==orbitals_index_gram );
    
    dm_->trmm('l', 'l', 't', 'n', 1., ls);
    dm_->trmm('r', 'l', 'n', 'n', 1., ls);
    
    uniform_occ_=false;
    occ_uptodate_=false;
    stripped_    =true;
}

void DensityMatrix::dressUpS(const dist_matrix::DistMatrix<DISTMATDTYPE>& ls,
              const int new_orbitals_index)
{
    assert( stripped_ );
    
    ls.trtrs('l', 't', 'n', *dm_);
    work_->transpose(*dm_);
    *dm_=*work_;
    ls.trtrs('l', 't', 'n', *dm_);

    orbitals_index_=new_orbitals_index;
    occ_uptodate_=false;
    uniform_occ_=false;    
    stripped_    =false;
}

// dm -> S*dm*S
//void DensityMatrix::surroundByS(const dist_matrix::DistMatrix<DISTMATDTYPE>& matS,
//                                const int new_orbitals_index)
//{
//    assert( !stripped_ );
//
//    work_->symm('r', 'l', 1., matS, *dm_, 0.);
//    dm_->symm('l', 'l', 1., matS, *work_, 0.);
//    
//    orbitals_index_=new_orbitals_index;
//    occ_uptodate_=false;
//    stripped_    =false;
//}
