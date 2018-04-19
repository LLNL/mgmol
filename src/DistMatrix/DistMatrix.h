// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

////////////////////////////////////////////////////////////////////////////////
//
// DistMatrix.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: DistMatrix.h 22 2012-12-05 21:55:18Z jeanluc $

#ifndef DISTMATRIX_H
#define DISTMATRIX_H

#include "MGmol_blas1.h"
#include "Timer.h"
#include "mputils.h"

#include <cassert>
#include <complex>
#include <vector>
#include <string.h>

#ifdef ADD_
#define numroc     numroc_
#define pdswap     pdswap_
#define psswap     psswap_
#define pdlange    pdlange_
#define pslange    pslange_
#define pdlacp3    pdlacp3_
#define pslacp3    pslacp3_
#define pdtran     pdtran_
#define pstran     pstran_
#define pdsymm     pdsymm_
#define pssymm     pssymm_
#define pdgemm     pdgemm_
#define pdgemv     pdgemv_
#define pdsymv     pdsymv_
#define psgemm     psgemm_
#define psgemv     psgemv_
#define pssymv     pssymv_
#define pdsyrk     pdsyrk_
#define pssyrk     pssyrk_
#define pdgemr2d   pdgemr2d_
#define psgemr2d   psgemr2d_
#define pigemr2d   pigemr2d_
#define pdtrmm     pdtrmm_
#define pstrmm     pstrmm_
#define pdtrsm     pdtrsm_
#define pstrsm     pstrsm_
#define pdtrtrs    pdtrtrs_
#define pstrtrs    pstrtrs_
#define pdpotrf    pdpotrf_
#define pspotrf    pspotrf_
#define pdpotrs    pdpotrs_
#define pspotrs    pspotrs_
#define pdgetrf    pdgetrf_
#define psgetrf    psgetrf_
#define pdgetrs    pdgetrs_
#define psgetrs    psgetrs_
#define pdpotri    pdpotri_
#define pspotri    pspotri_
#define pdtrtri    pdtrtri_
#define pstrtri    pstrtri_
#define pdpocon    pdpocon_
#define pspocon    pspocon_
#define pdsygst    pdsygst_
#define pssygst    pssygst_
#define pdsyev     pdsyev_
#define pssyev     pssyev_
#define pdlatra    pdlatra_
#define pslatra    pslatra_
#define pdlaset    pdlaset_
#define pslaset    pslaset_
#define pdgesvd    pdgesvd_
#define psgesvd    psgesvd_

#define dlange     dlange_
#define slange     slange_
#define dscal      dscal_
#define dcopy      dcopy_
#define ddot       ddot_
#define daxpy      daxpy_
#define dsymm      dsymm_
#define ssymm      ssymm_
#define dgemm      dgemm_
#define dgemv      dgemv_
#define dsymv      dsymv_
#define sgemm      sgemm_
#define sgemv      sgemv_
#define ssymv      ssymv_
#define dsyrk      dsyrk_
#define dtrmm      dtrmm_
#define strmm      strmm_
#define dtrsm      dtrsm_
#define strsm      strsm_
#define dtrtrs     dtrtrs_
#define strtrs     strtrs_
#define dpotrf     dpotrf_
#define spotrf     spotrf_
#define dpotrs     dpotrs_
#define spotrs     spotrs_
#define dgetrf     dgetrf_
#define sgetrf     sgetrf_
#define dgetrs     dgetrs_
#define sgetrs     sgetrs_
#define dpotri     dpotri_
#define spotri     spotri_
#define dtrtri     dtrtri_
#define strtri     strtri_
#define dpocon     dpocon_
#define spocon     spocon_
#define dsygst     dsygst_
#define ssygst     ssygst_
#define dsyev      dsyev_
#define ssyev      ssyev_
#define dgesvd     dgesvd_
#define sgesvd     sgesvd_
#endif

typedef const char* const  Pchar;
typedef const int* const  Pint;
typedef const double* const  Pdouble;
typedef const float* const  Pfloat;

extern "C"
{
  // BLAS1
  void dscal(Pint, Pdouble, double*, Pint);
  double ddot(Pint, Pdouble, Pint, Pdouble, Pint);
  void sscal(Pint, Pfloat, float*, Pint);
  float sdot(Pint, Pfloat, Pint, Pfloat, Pint);
}

#if USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

typedef double DISTMATDTYPE;
//typedef float DISTMATDTYPE;

namespace dist_matrix{

class BlacsContext;


template <class T>
class DistMatrix
{
private:
    static Timer   matgather_tm_;
    static Timer   pdgemr2d_tm_;
    static Timer   potri_tm_;
    static Timer   potrf_tm_;
    static Timer   trtri_tm_;
    
    static int distmatrix_def_block_size_;
    static BlacsContext* default_bc_;
    
    std::string object_name_;
  
    const MPI_Comm    comm_global_;
    const BlacsContext& bc_;
    int ictxt_;
    int lld_;  // leading dimension of local matrix
    int m_, n_;   // size of global matrix
    int mb_, nb_; // size of blocks
    int mloc_, nloc_;   // size of local array
    int nprow_, npcol_; // number of process rows and cols in Blacs context
    int myrow_, mycol_; // position of my process in the process grid
    int mblocks_, nblocks_; // number of local blocks
    int desc_[9];
    bool active_;
    bool m_incomplete_, n_incomplete_; // this process has an incomplete block

    // local array
    T* val_;
    int size_;
    
    void setDiagonalValues(const T* const dmat);

public:
  
    std::string name()const{ return object_name_; }
    
    static void setBlockSize(const int nb)
    {
        distmatrix_def_block_size_=nb;
    }
    static int getBlockSize()
    {
        return distmatrix_def_block_size_;
    }
    
    static void setDefaultBlacsContext(BlacsContext* bc)
    {
        default_bc_=bc;
    }
    
    static Timer matgather_tm()
    {
        return matgather_tm_;
    }  
    static Timer pdgemr2d_tm()
    {
        return pdgemr2d_tm_;
    }  

    static Timer potri_tm()
    {
        return potri_tm_;
    }  
    static Timer trtri_tm()
    {
        return trtri_tm_;
    }  

    static Timer potrf_tm()
    {
        return potrf_tm_;
    }  

    int nprow()const{return nprow_; }
    int npcol()const{return npcol_; }
    const T val(const int i) const{ return val_[i]; }
    void setval(const int i, const T val)
    {
      assert( i<size_ );
      val_[i]=val; 
    }
    void setval(const int l, const int m, const int x, const int y, 
                const T val)
    { 
      val_[l*mb_+x+(m*nb_+y)*mloc_] = val; 
    }
    void addval(const int index, 
                const T val)
    { 
      assert( index<size_ );
      val_[index] += val; 
    }

    void addval(const int l, const int m, const int x, const int y, 
                const T val)
    { 
      //cout<<"Add "<<val<<" to element "<<l<<","<<m<<","<<x<<","<<y
      //    <<" on PE "<<myrow_<<","<<mycol_<<endl;
      assert( l*mb_+x+(m*nb_+y)*mloc_<size_ );
      val_[l*mb_+x+(m*nb_+y)*mloc_] += val; 
    }

    //input: global index (i,j)
    void addval(const int i, const int j, 
                const T val)
    { 
      const int ib=(i/(nprow_*mb_));
      const int x=i % mb_;

      const int jb=(j/(npcol_*nb_));
      const int y=j % nb_;

      assert( ib*mb_+x+(jb*nb_+y)*mloc_<size_ );
      val_[ib*mb_+x+(jb*nb_+y)*mloc_] += val; 
    }
    
    int ictxt(void) const { return ictxt_; }
    int lld(void) const { return lld_; }
    int m(void) const { return m_; } // size of global matrix
    int n(void) const { return n_; } // size of global matrix
    int mb(void) const { return mb_; } // size of blocks
    int nb(void) const { return nb_; } // size of blocks
    int mloc(void) const { return mloc_; } // size of local array
    int nloc(void) const { return nloc_; } // size of local array
    int localsize(void) const { return mloc_*nloc_; } // local size of val
    size_t memsize(void) const { return m_ * n_ * sizeof(T); }
    size_t localmemsize(void) const { return mloc_ * nloc_ * sizeof(T); }
    const int* desc(void) const{ return &desc_[0]; }
    
    // local block size of block (l,m)
    int mbs(const int l) const
    {
      // if l is the last block and this process holds an incomplete
      // block, then the size is m_%mb_. Otherwise, it is mb_.
      return ( (l==mblocks_-1) && ( m_incomplete_ ) ) ? m_%mb_ : mb_;
    }
    int nbs(const int m) const
    {
      // if m is the last block and this process holds an incomplete
      // block, then the size is n_%nb_. Otherwise, it is nb_.
      return ( (m==nblocks_-1) && ( n_incomplete_ ) ) ? n_%nb_ : nb_;
    }
    
    // number of local blocks
    int mblocks(void) const { return mblocks_; }     
    int nblocks(void) const { return nblocks_; }
     
    int myrow(void) const { return myrow_; }     
    int mycol(void) const { return mycol_; }
     
    // functions to compute local indices from global indices

    // index of blocks: element (i,j) is at position (x,y)
    // in local block (l,m) of process (pr,pc)    
    int ib(const int i) const { 
      assert( i<m_ );
      assert( i>=0 );
      return i/(nprow_*mb_); 
    }
    int x(const int i) const {
      assert( i<m_ );
      assert( i>=0 );
      return i % mb_;
    }
    int pr(const int i) const {
      assert( i<m_ );
      assert( i>=0 );
      return (i/mb_) % nprow_;
    }
    
    int jb(const int j) const {
      assert( j<n_ );
      assert( j>=0 );
      return j/(npcol_*nb_);
    }
    int y(const int j) const {
      assert( j<n_ );
      assert( j>=0 );
      return j % nb_;
    }
    int pc(const int j) const {
      assert( j<n_ );
      assert( j>=0 );
      return (j/nb_) % npcol_;
    }
    
    // global indices:
    // (i,j) is the global index of element (x,y) of block (l,m)
    int i(const int l, const int x) const
    {
        return (l * nprow_ + myrow_) * mb_ + x;
    }     
    int j(const int m, const int y) const
    {
        return (m * npcol_ + mycol_) * nb_ + y;
    }
    
    int indxl2grow(const int indxloc)const;
    int indxl2gcol(const int indxloc)const;

    // active flag: the matrix has elements on this process
    bool active(void) const { return active_; }
    
    DistMatrix<T>(const std::string& name);
        
    DistMatrix<T>(const std::string& name, const BlacsContext& bc);
        
    // Construct a DistMatrix of dimensions m,n
    DistMatrix<T>(const std::string& name, const BlacsContext& bc, 
                  const int m, const int n);
                  
    DistMatrix<T>(const std::string& name,
                  const int m, const int n);
                  
    // Construct a DistMatrix of dimensions m,n
    DistMatrix<T>(const std::string& name, const BlacsContext& bc, 
                  const int m, const int n,
                  const int mb, 
                  const int nb);
                  
    DistMatrix<T>(const std::string& name, 
                  const int m, const int n,
                  const int mb, 
                  const int nb);
                  
    // copy constructor
    DistMatrix<T>(const DistMatrix<T>& a);

    // Construct a DistMatrix from a duplicated matrix a, with leading dim lda
    DistMatrix<T>(const std::string& name,
                  const BlacsContext& bc, 
                  const int lda, const T* const a, 
		  const int m, const int n,
                  const int mb, const int nb);
                  
    DistMatrix<T>(const std::string& name,
                  const BlacsContext& bc, 
                  const int lda, const T* const a, 
		  const int m, const int n);
                  
    // Construct a diagonal DistMatrix from a vector dmat of diagonal elements
    DistMatrix<T>(const std::string& name,
                  const BlacsContext&,
    		  const T * const dmat,
    		  const int m, const int n,
    		  const int mb, const int nb);
                  
    // Construct a diagonal DistMatrix from a vector dmat of diagonal elements
    DistMatrix<T>(const std::string& name,
                  const BlacsContext&,
    		  const T * const dmat,
    		  const int m, const int n);
                  
    DistMatrix<T>& operator=(const DistMatrix<T>& a);
    DistMatrix<T>& assign(const DistMatrix<T>&, const int, const int);
    DistMatrix<T>& operator+=(const DistMatrix<T>& a) { axpy(1.0,a); return *this; }
    DistMatrix<T>& operator-=(const DistMatrix<T>& a) { axpy(-1.0,a); return *this; }

    ~DistMatrix()
    {
      if( active_ )
         delete[] val_;
    }

    void identity(void);
    void setVal(const char, const T);
    void matgather(T *a, const int lda)const;

    void resize(const int m, const int n, 
                const int mb, const int nb);    
    void init(const T* const a, const int lda);
    void initTest();

    void clear(void);
    
    double dot(const DistMatrix<T> &a)const;
    void axpy(const double alpha, const DistMatrix<T> &a);
    void scal(const double alpha);
    double nrm2()const;
    T asum(DistMatrix<T> &a)const;
    T amax(DistMatrix<T> &a)const;
    
    double trace(void)const;
    void setDiagonal(const std::vector<T>& diag_values);

    // matrix * matrix
    // this = alpha*op(A)*op(B)+beta*this
    void gemm(const char transa, const char transb,
         const T alpha, const DistMatrix<T>& a,
         const DistMatrix<T>& b, const T beta);
    
    // symmetric_matrix * matrix
    // this = alpha * A * B + beta * this
    void symm(const char side, const char uplo,
         const T alpha, const DistMatrix<T>& a,
         const DistMatrix<T>& b, const T beta);
    
    // matrix * vector
    // this = beta*this + alpha * A * x
    void gemv(const char trans, const T alpha, const DistMatrix<T>& a,
         const DistMatrix<T>& x, const T beta);
    // symmetric version of matvec
    void symv(const char uplo, const T alpha, const DistMatrix<T>& a,
         const DistMatrix<T>& x, const T beta);

    // symmetric rank k update
    // this = beta * this + alpha * A * A^T  (trans=='n')
    // this = beta * this + alpha * A^T * A  (trans=='t')
    void syrk(char uplo, char trans,
         T alpha, DistMatrix<T>& a, T beta);

    // matrix transpose
    // this = alpha * transpose(A) + beta * this
    void transpose(const T alpha,const  DistMatrix<T>& a, const T beta);
    void transpose(const DistMatrix<T>& a);

    void trset(const char);

    void trmm(const char, const char, const char, const char, 
              const T, const DistMatrix<T>&);
    void trsm(const char, const char,const  char, const char, 
              const T, const DistMatrix<T>&);
    void trtrs(const char, const char, const char, DistMatrix<T>&)const;

    // Cholesky decomposition of a symmetric matrix
    int potrf(char uplo);
    // Inverse of a symmetric matrix from Cholesky factor
    int potri(char uplo);
    void potrs(char, DistMatrix<T>&);
    
    int trtri(char uplo, char diag);

    void getrf(std::vector<int>&);
    void getrs(char, DistMatrix<T>&, std::vector<int>&);
    double norm(char);
    double pocon(char,T);
    void sygst(int, char, const DistMatrix<T>&);
    void syev(char, char, std::vector<T>&, DistMatrix<T>&);
    void sygv(char, char, const char, DistMatrix<T>&, std::vector<T>&, DistMatrix<T>&);
    void gesvd(char jobu, char jobvt, std::vector<T>& s,
               DistMatrix<T>& u,
               DistMatrix<T>& v);
    
    void getsub(const DistMatrix<T> &a, int m, int n, int ia, int ja);
    
    int iamax(const int j, T& val);

    void print(std::ostream& os)const;
    void print(std::ostream& os, const int, const int, const int, const int) const;
    void printMM(std::ostream& os)const;

    // get value for global index (i,j)
    // assuming we are on the right processor to get it!
    const T getGlobalVal(const int i, const int j,
                         const bool onpe=false)const;
    
    // set value for global index (i,j)
    void setGlobalVal(const int i, const int j, const T val,
                      const bool onpe=false);
    
    void assign(const T* const v, const int n)
    {
      if( active_ ){
        assert( n<=size_ );
        memcpy(val_, v, n*sizeof(T) );
      }
    }
    void assignColumn(const T* const v, const int i)
    {
      if( active_ ){
        assert( i<=nloc_ );
        if(size_>0)memcpy(&val_[i*mloc_], v, mloc_*sizeof(T) );
      }
    }
    void swapColumns(const int j1, const int j2);

    void axpyColumn(const int icol, const double alpha, const DistMatrix<T>& x);
    void scalColumn(const int icol, const double alpha);
    
    void copyDataToVector(std::vector<T>& v)const
    {
      assert( (int)v.size()>=size_);
      memcpy(&v[0],&val_[0],mloc_*nloc_*sizeof(T));
    }
    double dotColumns(const int i,const int j)const
    {
//      int ione=1;
      return MPdot(mloc_,&val_[mloc_*i],&val_[mloc_*j]);
//      return ddot(&mloc_,&val_[mloc_*i],&ione,&val_[mloc_*j],&ione);
    }
    double dotColumns(const int i,
                      const DistMatrix<T> &a, const int j)const
    {
//      int ione=1;
      return MPdot(mloc_,&val_[mloc_*i],&a.val_[mloc_*j]);
//      return ddot(&mloc_,&val_[mloc_*i],&ione,&a.val_[mloc_*j],&ione);
    }
    double sumProdElements(const DistMatrix<T>& a)const;
    void getDiagonalValues(T* const dmat)const;
    
    MPI_Comm comm_global()const;
};
template <class T> std::ostream& operator << ( std::ostream& os, DistMatrix<T>& a );

} // namespace

#endif
