// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
 * Main header file for data distribution
*/

#ifndef _DATADISTRIBUTION_H_
#define _DATADISTRIBUTION_H_

#include "VariableSizeMatrix.h"
#include "PEenv.h"
#include "MPIdata.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>

class DataDistribution
{
  static Timer gathersizes_tm_;
  static Timer reducesizes_tm_;
  static Timer distribute_local_row_tm_;
  static Timer distribute_local_data_tm_;
  static Timer send_recv_tm_;  
  static Timer send_recv_rows_tm_;  
  static Timer update_row_tm_;  
  static Timer augment_local_data_tm_;
  static Timer update_local_rows_tm_;
  static Timer initLocalRow_tm_;
  static Timer send_recv_ovlp_tm_;
  static Timer send_recv_rows_ovlp_tm_;
  static Timer send_recv_ovlp_wait_tm_;
  static Timer send_recv_rows_ovlp_wait_tm_;  
  
  static int max_matsize_;
  static int max_nnz_;
  
  static short count_computeMaxDataSize_;
  static short maxcount_computeMaxDataSize_;
  
  std::string name_;
  const double spread_radius_;				// Spreading radius for data distribution
  const pb::PEenv& mypeenv_;

  int lsize_;				// Initial size of local matrix
  int aug_size_;				// Augmented size of the local matrix after data distribution
  double domain_[3];			// Problem domain dimensions
  int nproc_xyz_[3];			// Number of processors in each direction;
  double proc_width_[3];				// Processor radius in each direction  
  MPI_Comm cart_comm_;				// MPI cartesian communicator for data distribution
  MPI_Request request_[4];	// MPI request handle
  int sbuf_[2]; 		// work buffer for data transfer (send)
  int rbuf_[2];		// work buffer for data transfer (recv)
  int lstep_[3];			// number of steps to the left for each dimension
  int rstep_[3];			// number of steps to the right for each dimension
  int data_pos_[6]; 		// positions of datatypes within buffer - 2 entries per direction	
  int loc_data_sz_[3];		// size of initial local data that is packed for data distribution

/* pointers to communication recv buffer entries */
  /* position of number of rows of matrix in recv buffer */
  int* rbuf_nrows_ptr_;
  /* start position of double datatypes in buffer */
  int* rbuf_start_double_pos_ptr_;
  /* start position of matrix nonzero row count info */
  int* rbuf_nnzrow_ptr_;    
  /* start position of local variables global indexes array info */
   int* rbuf_lvars_ptr_;
  /* start position of column index array info */    
  int* rbuf_pj_ptr_; 
  /* start position of matrix (double) coefficient data */    
  double* rbuf_pa_ptr_;
  /* actual size of data in recv buffer (in bytes or sizeof char) */
  int rbuf_data_size_;

/* Compute number of steps in each direction */
  void computeNumSteps() ; 
/* compute number of steps given the max number of steps in each direction */
  void computeNumSteps(const int max_steps[3]) ; 
//  template <class T>
  void gatherDataSizes(const short dir, const VariableSizeMatrix<sparserow>& lmat, int *maxsize, int *nzmax); /* compute max data sizes from neighbors*/
  /* compute starting positions for packing local data */
//  template <class T>
  void computePackedDataPositions(const VariableSizeMatrix<sparserow>& lmat, int *pos)const
  {
      const int lsize = lmat.n();
      const int nnzmat = lmat.nnzmat();
   
      /* now compute start positions for packing local data */
      pos[0] = 0; /* position of integer variables */ 
      const int offset = 2*(lsize+1)*sizeof(int) + (nnzmat+1)*sizeof(int);
      const int padding = (offset%sizeof(double));
      pos[1] = offset + padding;   /* position of double variables - adjusted for alignment */
      
      return;
  }
  /* get maximum buffer size for data transfer */
  int getPackedBufferSize(const int maxsize, const int nzmax)const
  {
      /* compute buffer size for packing data - include alignment padding for double, if any */
      const int offset = 2*(maxsize + 1)*sizeof(int) + (nzmax+1)*sizeof(int);
      const int padding = (offset%sizeof(double));
      const int bsiz = offset + nzmax*sizeof(double) + padding;
    
      return bsiz;
  }

  /* setup persistent requests */
  void setupPersistentRequests(const short dir)
  {
    /* initialize persisted requests */
    request_[0] = MPI_REQUEST_NULL;
    request_[1] = MPI_REQUEST_NULL;    
    request_[2] = MPI_REQUEST_NULL;
    request_[3] = MPI_REQUEST_NULL; 
    
     int source, dest;
     /* Get source and destination ID to send and recv data - left direction */
     short disp = -1;
     MPI_Cart_shift(cart_comm_, dir, disp, &source, &dest);
     /* setup request for left direction */
     MPI_Recv_init(&rbuf_[0], 2, MPI_INT, source, 0, cart_comm_, &request_[0]);
     MPI_Rsend_init(&sbuf_[0], 2, MPI_INT, dest, 0, cart_comm_, &request_[1]); 

     /* setup request for right direction - just reverse source and dest (no need for mpi_cart_shift()) */
     MPI_Recv_init(&rbuf_[0], 2, MPI_INT, dest, 1, cart_comm_, &request_[2]);
     MPI_Rsend_init(&sbuf_[0], 2, MPI_INT, source, 1, cart_comm_, &request_[3]);      
     
     return;
  }
  /* delete persistent requests */
  void deletePersistentRequests()
  {
       MPI_Request_free(&request_[0]);
       MPI_Request_free(&request_[1]);
       MPI_Request_free(&request_[2]);       
       MPI_Request_free(&request_[3]);  
       
       return;
  }
  
  /* Set data pointer positions on recv buffer */
  void setPointersToRecvData(const char *rbuf)
  {
     assert(rbuf != NULL);
     
     rbuf_nrows_ptr_ = (int *)rbuf;
     /* check if matrix is nonzero */
     if(*rbuf_nrows_ptr_ != 0)
     {
        rbuf_start_double_pos_ptr_ = rbuf_nrows_ptr_ + 1;
        rbuf_nnzrow_ptr_ = rbuf_start_double_pos_ptr_ + 1;
        rbuf_lvars_ptr_ = rbuf_nnzrow_ptr_ + (*rbuf_nrows_ptr_ + 1);
        rbuf_pj_ptr_ = rbuf_lvars_ptr_ + (*rbuf_nrows_ptr_);
        rbuf_pa_ptr_ = (double *)(rbuf + *rbuf_start_double_pos_ptr_);
        
        /* compute size of data in buffer */
        rbuf_data_size_ = *rbuf_start_double_pos_ptr_ + rbuf_nnzrow_ptr_[*rbuf_nrows_ptr_]*sizeof(double);
     }
     else
     {
        /* The size is equal to the size of an int if nrows == 0. */
        rbuf_data_size_ = sizeof(int);     
     }
  }

  /* Reset data pointer positions on recv buffer to NULL*/
  void resetPointersToRecvDataToNULL()
  {
     rbuf_nrows_ptr_ = NULL;
     rbuf_start_double_pos_ptr_ = NULL;
     rbuf_nnzrow_ptr_ = NULL;
     rbuf_lvars_ptr_ = NULL;
     rbuf_pj_ptr_ = NULL;
     rbuf_pa_ptr_ = NULL;
        
     rbuf_data_size_ = 0;
  }
  
  /* Perform data distribution of local data */  
  template <class T>
  void distributeLocalData(const int nsteps, const int dir, const int disp, const int bsiz, const int *pos, 
                           VariableSizeMatrix<sparserow>& lmat, VariableSizeMatrix<T>& aug_mat, const bool append, 
                           const bool bcflag);  
  template <class T>
  void distributeLocalDataWithCommOvlp(const int nsteps, const int dir, const int disp, const int bsiz, const int *pos, 
                           VariableSizeMatrix<sparserow>& lmat, VariableSizeMatrix<T>& aug_mat, const bool append, 
                           const bool bcflag);
  template <class T>  
  void distributeLocalRows(const int nsteps, const int dir, const int disp, const int bsiz, const int *pos, 
                           VariableSizeMatrix<sparserow>& lmat, VariableSizeMatrix<T>& aug_mat, const bool append);
  template <class T>  
  void distributeLocalRowsWithCommOvlp(const int nsteps, const int dir, const int disp, const int bsiz, const int *pos, 
                           VariableSizeMatrix<sparserow>& lmat, VariableSizeMatrix<T>& aug_mat, const bool append);                           
  template<class T>
  void mergeDataFromNeighborToLocalData(VariableSizeMatrix<T>& aug_mat, const char *rbuf, const bool append);
  template<class T>
  void updateExistingLocalDataEntriesWithRecvBuf(VariableSizeMatrix<T>& amat, const char *rbuf);
  template<class T>
  void copyRowsFromRecvBuf(VariableSizeMatrix<T>& amat, const char *rbuf, const bool append);  
  /* Pack local data into buffer */
  void packLocalData(VariableSizeMatrix<sparserow>& lmat, const int *pos, char *buf);  
  /* get recv buffer data size */
  int getRecvBufDataSize(){return rbuf_data_size_;}
  void reduceDataSizes(const short dir, const VariableSizeMatrix<sparserow>& lmat, int *maxsize, int *nzmax);
  
  void computeMaxDataSize(const short dir, const VariableSizeMatrix<sparserow>& lmat, int *maxsize, int *nzmax);
 
public:
    DataDistribution(const std::string name, const double s_radius, const pb::PEenv& myPEenv, const double domain[]);
    DataDistribution(const std::string name, const int max_steps[3], const pb::PEenv& myPEenv, const double domain[]);
    ~DataDistribution(){}
    
    static void enforceComputeMaxDataSize()
    {
        count_computeMaxDataSize_=maxcount_computeMaxDataSize_-1;
    }
  
    template <class T>
    void augmentLocalData(VariableSizeMatrix<T>& vsmat, const bool append, const bool bcflag=false);	// augment the local matrix
    
    template <class T>    
    void updateLocalRows(VariableSizeMatrix<T>& vsmat, const bool append=false);				// augment the local matrix
    
    static void printTimers(std::ostream& os); // print timers   
    
    void printStats()
    {
      if(onpe0)
      {
         std::cout<<"spread_radius = "<<spread_radius_<<" initial size = "<<lsize_<<" augmented size = "<<aug_size_<<std::endl;
         std::cout<<"x-direction."<<std::endl;
         std::cout<<"lstep = "<<lstep_[0]<<", rstep = "<<rstep_[0]<<std::endl; 
         std::cout<<"y-direction."<<std::endl;
         std::cout<<"lstep = "<<lstep_[1]<<", rstep = "<<rstep_[1]<<std::endl; 
         std::cout<<"z-direction."<<std::endl;
         std::cout<<"lstep = "<<lstep_[2]<<", rstep = "<<rstep_[2]<<std::endl;         
       }
       return;
    }
};

#endif  
