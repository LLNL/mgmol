// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: PEenv.cc,v 1.18 2009/08/31 16:22:51 jeanluc Exp $
#include "PEenv.h"
#include "tools.h"

#include <cassert>
#include <stdlib.h>
using namespace std;

namespace pb
{

PEenv::PEenv(MPI_Comm comm, ostream* os)
    : // default constructor
      comm_(comm),
      os_(os)
{
    // default values for 1 MPI job
    //(*os_)<<" default constructor for PEenv\n";

    mytask_      = 0;
    n_mpi_tasks_ = 1;
    for (int i = 0; i < 3; i++)
    {
        n_mpi_tasks_dir_[i] = 1;
        mytask_dir_[i]      = 0;
    }
    for (int i = 0; i < 6; i++)
        mpi_neighbors_[i] = 0;

    color_       = 0;
    comm_active_ = 0;
    cart_comm_   = 0;
    comm_x_      = 0;
    comm_y_      = 0;
    comm_z_      = 0;

    onpe0_ = true;
    for (int i = 0; i < 3; i++)
        other_tasks_dir_[i] = NULL;
#ifdef USE_MPI
    comm_active_ = MPI_COMM_SELF;
    int ltask;
    MPI_Comm_rank(comm_, &ltask);
    if (ltask > 0) onpe0_ = false;
        // cout<<"ltask="<<ltask<<", comm_="<<comm_<<endl;
#endif
}

// Constructor for a global grid nx by ny by nz
// Must be called from all the PEs simultaneously!
PEenv::PEenv(MPI_Comm comm, const int nx, const int ny, const int nz,
    const int bias, ostream* os)
    : comm_(comm), os_(os)
{
    assert(nx > 0);
    assert(ny > 0);
    assert(nz > 0);

    // default values for 1 MPI job
    // cout<<" default constructor for PEenv\n";
    mytask_      = 0;
    n_mpi_tasks_ = 1;
    for (int i = 0; i < 3; i++)
    {
        n_mpi_tasks_dir_[i] = 1;
        mytask_dir_[i]      = 0;
    }
    for (int i = 0; i < 6; i++)
        mpi_neighbors_[i] = 0;

    color_ = 0;
    onpe0_ = true;
#ifdef USE_MPI
    comm_active_ = comm_;

    int mpirc = MPI_Comm_size(comm_, &n_mpi_tasks_);
    mpirc     = MPI_Comm_rank(comm_, &mytask_);

    if (mytask_ > 0) onpe0_ = false;

    split_comm(nx, ny, nz, bias);

    int periods[3] = { 1, 1, 1 };
    int reorder    = 0;
    MPI_Cart_create(
        comm_active_, 3, n_mpi_tasks_dir_, periods, reorder, &cart_comm_);

    int mytask_cart_comm_;
    MPI_Comm_rank(cart_comm_, &mytask_cart_comm_);
    MPI_Comm_rank(comm_, &mytask_);
    assert(mytask_cart_comm_ == mytask_);

    MPI_Comm_size(comm_, &n_mpi_tasks_);
    // if(color_>0){
    //    n_mpi_tasks_=0;
    //    mytask_=-1;
    //    for(int i=0;i<3;i++){
    //        n_mpi_tasks_dir_[i]=0;
    //    }
    //}

    MPI_Barrier(comm_);

    task2xyz();
    setup_my_neighbors();

    MPI_Barrier(comm_);

    {
        // create a unique color based on y and z coordinate in processor grid
        int color = n_mpi_tasks_dir_[2] * mytask_dir_[1] + mytask_dir_[2];
        mpirc     = MPI_Comm_split(cart_comm_, color, mytask_dir_[0], &comm_x_);

        int task = -1;
        if (mytask_dir_[0] == 0)
        {
            task = mytask_;
        }
        MPI_Allreduce(&task, &root_x_, 1, MPI_INT, MPI_MAX, comm_x_);
    }
    {
        // create a unique color based on x and z coordinate in processor grid
        int color = n_mpi_tasks_dir_[2] * mytask_dir_[0] + mytask_dir_[2];
        mpirc     = MPI_Comm_split(cart_comm_, color, mytask_dir_[1], &comm_y_);
    }
    {
        // create a unique color based on x and y coordinate in processor grid
        int color = n_mpi_tasks_dir_[1] * mytask_dir_[0] + mytask_dir_[1];
        mpirc     = MPI_Comm_split(cart_comm_, color, mytask_dir_[2], &comm_z_);
    }

    MPI_Barrier(comm_);
#else
    (void)nx; // unused
    (void)ny; // unused
    (void)nz; // unused
    comm_active_ = 0;
#endif

    set_other_tasks_dir();
}

void PEenv::set_other_tasks_dir()
{
    for (int i = 0; i < 3; i++)
        other_tasks_dir_[i] = new int[n_mpi_tasks_];

#ifdef USE_MPI
    int* buffer = new int[3 * n_mpi_tasks_];
    MPI_Allgather(&mytask_dir_[0], 3, MPI_INT, buffer, 3, MPI_INT, comm());
    for (int j = 0; j < n_mpi_tasks_; j++)
    {
        for (int i = 0; i < 3; i++)
        {
            other_tasks_dir_[i][j] = buffer[3 * j + i];
        }
    }
    delete[] buffer;
#else
    for (int i = 0; i < 3; i++)
        other_tasks_dir_[i][0] = 0;
#endif
}

#if 0
PEenv::PEenv(MPI_Comm comm, 
             const int mytask, const int taskx, const int tasky, const int taskz,
             ostream* os):
    comm_(comm),
    os_(os)
{
    //cout<<" constructor for PEenv\n";

    mytask_=mytask;
    n_mpi_tasks_=taskx*tasky*taskz;
    for(int i=0;i<3;i++){
        n_mpi_tasks_dir_[i]=1;
        mytask_dir_[i]=0;
    }    
    for(int i=0;i<6;i++)
        mpi_neighbors_[i]=0;
    
    color_=0;
    onpe0_= (mytask_==0);
#ifdef USE_MPI
    comm_active_=comm;

    n_mpi_tasks_dir_[0]=taskx;
    n_mpi_tasks_dir_[1]=tasky;
    n_mpi_tasks_dir_[2]=taskz;

    task2xyz();
    setup_my_neighbors();
#else
    comm_active_=0;
    assert(n_mpi_tasks_==1);
    assert(mytask==0);
    assert(taskx==1);
    assert(tasky==1);
    assert(taskz==1);
#endif
}
#endif

PEenv::~PEenv()
{
#ifdef USE_MPI
    if (comm_active_ != comm_ && comm_active_ != MPI_COMM_SELF)
    {
        // cout<<"MPI_Comm_free: "<<comm_active_<<endl;
        assert(comm_active_ != NULL);
        int mpirc    = MPI_Comm_free(&comm_active_);
        comm_active_ = 0;
        if (mpirc != MPI_SUCCESS)
        {
            cerr << "MPI_Comm_free failed!" << endl;
            MPI_Abort(comm_, 2);
        }
    }
    if (cart_comm_ != 0) MPI_Comm_free(&cart_comm_);
    if (comm_x_ != 0) MPI_Comm_free(&comm_x_);
#endif
    for (int i = 0; i < 3; i++)
    {
        if (other_tasks_dir_[i] != NULL) delete[] other_tasks_dir_[i];
    }
}

void PEenv::barrier(void) const
{
#ifdef USE_MPI
    MPI_Barrier(comm_);
#endif
}

double PEenv::double_max_all(double x) const
{
#ifdef USE_MPI
    double out;

    int mpirc = MPI_Allreduce(&x, &out, 1, MPI_DOUBLE, MPI_MAX, comm_active_);
    if (mpirc != MPI_SUCCESS)
    {
        cerr << "real_max_all:MPI_Allreduce double max failed!!!" << endl;
        exit(1);
    }
    return out;
#else
    return x;
#endif
}

double PEenv::double_min_all(double x) const
{
#ifdef USE_MPI
    double out;

    int mpirc = MPI_Allreduce(&x, &out, 1, MPI_DOUBLE, MPI_MIN, comm_active_);
    if (mpirc != MPI_SUCCESS)
    {
        cerr << "real_max_all:MPI_Allreduce double max failed!!!" << endl;
        exit(1);
    }
    return out;
#else
    return x;
#endif
}

int PEenv::int_max_all(int x)
{
#ifdef USE_MPI
    int out;

    int mpirc = MPI_Allreduce(&x, &out, 1, MPI_INT, MPI_MAX, comm_active_);
    if (mpirc != MPI_SUCCESS)
    {
        cerr << "int_max_all: MPI_Allreduce int max failed!!!" << endl;
        globalExit(2);
    }
    return out;
#else
    return x;
#endif
}

double PEenv::double_sum_all(double x) const
{
#ifdef USE_MPI
    double out;

    int mpirc = MPI_Allreduce(&x, &out, 1, MPI_DOUBLE, MPI_SUM, comm_active_);
    if (mpirc != MPI_SUCCESS)
    {
        cerr << "double_sum_all: MPI_allreduce double sum failed!!!" << endl;
        globalExit(2);
    }
    return out;

#else

    return x;

#endif
}

void PEenv::task2xyz()
{
    if (color_ == 0)
    {
        assert(n_mpi_tasks_dir_[0] > 0);
        assert(n_mpi_tasks_dir_[1] > 0);
        assert(n_mpi_tasks_dir_[2] > 0);
#if 1
        int rc = MPI_Cart_coords(cart_comm_, mytask_, 3, mytask_dir_);
        if (rc != MPI_SUCCESS)
        {
            cerr << " error in MPI_Cart_coords()!!!" << endl;
            MPI_Abort(comm_, 1);
        }
#else
        mytask_dir_[2] = mytask_ % n_mpi_tasks_dir_[2];

        mytask_dir_[0] = mytask_ / n_mpi_tasks_dir_[2];
        mytask_dir_[1] = mytask_dir_[0] % n_mpi_tasks_dir_[1];
        mytask_dir_[0] /= n_mpi_tasks_dir_[1];

        if (mytask_dir_[0] >= n_mpi_tasks_dir_[0])
            mytask_dir_[0] -= n_mpi_tasks_dir_[0];
        if (mytask_dir_[0] >= n_mpi_tasks_dir_[0])
            mytask_dir_[0] -= n_mpi_tasks_dir_[0];
#endif
        assert(my_mpi(0) >= 0);
        assert(my_mpi(1) >= 0);
        assert(my_mpi(2) >= 0);

        // cout<<" task "<<mytask_<<" ->
        // ("<<my_mpi(0)<<","<<my_mpi(1)<<","<<my_mpi(2)<<")\n";
    }
}

void PEenv::setup_my_neighbors()
{
    assert(cart_comm_ != 0);

    if (color_ == 0)
    {
#if 1
        MPI_Cart_shift(
            cart_comm_, 0, 1, &mpi_neighbors_[WEST], &mpi_neighbors_[EAST]);
        MPI_Cart_shift(
            cart_comm_, 1, 1, &mpi_neighbors_[SOUTH], &mpi_neighbors_[NORTH]);
        MPI_Cart_shift(
            cart_comm_, 2, 1, &mpi_neighbors_[DOWN], &mpi_neighbors_[UP]);
#else
        mpi_neighbors_[EAST]  = xyz2task(my_mpi(0) + 1, my_mpi(1), my_mpi(2));
        mpi_neighbors_[WEST]  = xyz2task(my_mpi(0) - 1, my_mpi(1), my_mpi(2));
        mpi_neighbors_[NORTH] = xyz2task(my_mpi(0), my_mpi(1) + 1, my_mpi(2));
        mpi_neighbors_[SOUTH] = xyz2task(my_mpi(0), my_mpi(1) - 1, my_mpi(2));
        mpi_neighbors_[UP]    = xyz2task(my_mpi(0), my_mpi(1), my_mpi(2) + 1);
        mpi_neighbors_[DOWN]  = xyz2task(my_mpi(0), my_mpi(1), my_mpi(2) - 1);
#endif
        for (int i = 0; i < 6; i++)
        {
            if (mpi_neighbors(i) >= n_mpi_tasks())
            {
                cout << "mpi_neighbors(" << i << ")=" << mpi_neighbors(i)
                     << endl;
                cout << "n_mpi_tasks()=" << n_mpi_tasks() << endl;
            }
            assert(mpi_neighbors(i) >= 0);
            assert(mpi_neighbors(i) < n_mpi_tasks());
        }
    }
}

#ifdef USE_MPI

int PEenv::geom(const int nx, const int ny, const int nz, const int bias)
{
    int n[4] = { nx, ny, nz, n_mpi_tasks_ };

    const int pmax        = 9;
    const int coeff[pmax] = { 2, 3, 5, 7, 11, 13, 17, 19, 23 };

    int fac[4][pmax + 2];

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < pmax + 1; j++)
            fac[i][j] = 0;

    if (onpe0_ && os_ != 0)
        (*os_) << "Factorization of mesh dimensions and number of cpus:"
               << endl;
    for (int i = 0; i < 4; i++)
    {
        const int ni = n[i];
        for (int p = 0; p < pmax; p++)
        {
            const int nfactor = coeff[p];
            for (int k = 0; k < 20; k++)
            {
                if ((n[i] % nfactor) == 0)
                {
                    n[i] /= nfactor;
                    fac[i][p]++;
                }
                else
                {
                    break;
                }
            }
        }

        fac[i][pmax] = n[i];
        if (onpe0_ && os_ != 0)
            (*os_) << ni << "="
                   << " 2^" << fac[i][0] << "*3^" << fac[i][1] << "*5^"
                   << fac[i][2] << "*7^" << fac[i][3] << "*11^" << fac[i][4]
                   << "*13^" << fac[i][5] << "*17^" << fac[i][6] << "*19^"
                   << fac[i][7] << "*23^" << fac[i][8] << "*" << fac[i][9]
                   << endl;
    }

    for (int i = 0; i < 3; i++)
    {
        if (fac[i][0] > 1)
        {
            fac[i][0] -= 2;
        }
        else
        {
            if (onpe0_ && os_ != 0)
            {
                (*os_) << "Direction " << i << ": ";
                (*os_) << "Poisson Solver Requires grid size to be divisible "
                          "by 4!!!"
                       << endl;
            }
            return 0;
        }
    }

    // reset n[]
    n[0] = nx;
    n[1] = ny;
    n[2] = nz, n[3] = n_mpi_tasks_;
#ifndef NDEBUG
    if (onpe0_ && os_ != 0)
        (*os_) << "n=" << n[0] << ", " << n[1] << ", " << n[2] << ", " << n[3]
               << endl;
#endif

    bool div = true;
#if 0 // if we want to divide by factor 2 first
    // pick the largest factor 2
    while( fac[3][0]>0 && div){
        
        div=false;
        
        int imax0=2;
        int imax1=1;
        int imax2=0;
        if( (fac[2][0]>=fac[1][0])&&(fac[1][0]>=fac[0][0]) )
        {
            imax0=2; imax1=1; imax2=0;
        }
        else if( (fac[1][0]>=fac[2][0])&&(fac[2][0]>=fac[0][0]) )
        {
            imax0=1; imax1=2; imax2=0;
        }
        else if( (fac[2][0]>=fac[0][0])&&(fac[0][0]>=fac[1][0]) )
        {
            imax0=2; imax1=0; imax2=1;
        }
        else if( (fac[1][0]>=fac[0][0])&&(fac[0][0]>=fac[2][0]) )
        {
            imax0=1; imax1=0; imax2=2;
        }
        else if( (fac[0][0]>=fac[2][0])&&(fac[2][0]>=fac[1][0]) )
        {
            imax0=0; imax1=2; imax2=1;
        }
        else if( (fac[0][0]>=fac[1][0])&&(fac[1][0]>=fac[2][0]) )
        {
            imax0=0; imax1=1; imax2=2;
        }
        
        // Try first with the largest factor 2
        // Do imax0 first...
        if( fac[imax0][0]>0 ){
#ifndef NDEBUG
            if( onpe0_ && os_!=0 )
                (*os_)<<"Divide domain along direction "<<imax0<<endl;
#endif 
            n_mpi_tasks_dir_[imax0]*=2;
            fac[imax0][0]--;
            n[imax0]/=2;
            fac[3][0]--;
            n[3]/=2;
            div=true;
        }
       
        if(!div)
        if( fac[imax1][0]>0 ){
#ifndef NDEBUG
            if( onpe0_ && os_!=0 )
                (*os_)<<"Divide domain along direction "<<imax1<<endl;
#endif 
            n_mpi_tasks_dir_[imax1]*=2;
            fac[imax1][0]--;
            n[imax1]/=2;
            fac[3][0]--;
            n[3]/=2;
            div=true;
        }

        if(!div)
        if( fac[imax2][0]>0 ){
#ifndef NDEBUG
        if( onpe0_ && os_!=0 )
            (*os_)<<"Divide domain along direction "<<imax2<<endl;
#endif
            n_mpi_tasks_dir_[imax2]*=2;
            fac[imax2][0]--;
            n[imax2]/=2;
            fac[3][0]--;
            n[3]/=2;
            div=true;
        }
    }
#endif

    //    const int bias=4; // to change weight in dimension

    div = true;
    // Try to  divide the grid among PEs
    while (div)
    {
        div = false;

        // pick the largest grid dimension (with bias)
        int imax0 = 1;
        int imax1 = 2;
        int imax2 = 0;
        if ((n[2] >= n[1]) && (bias * n[1] >= n[0]))
        {
            imax0 = 2;
            imax1 = 1;
            imax2 = 0;
        }
        else if ((n[1] >= n[2]) && (bias * n[2] >= n[0]))
        {
            imax0 = 1;
            imax1 = 2;
            imax2 = 0;
        }
        else if ((bias * n[2] >= n[0]) && (n[0] >= bias * n[1]))
        {
            imax0 = 2;
            imax1 = 0;
            imax2 = 1;
        }
        else if ((bias * n[1] >= n[0]) && (n[0] >= bias * n[2]))
        {
            imax0 = 1;
            imax1 = 0;
            imax2 = 2;
        }
        else if ((n[0] >= bias * n[2]) && (n[2] >= n[1]))
        {
            imax0 = 0;
            imax1 = 2;
            imax2 = 1;
        }
        else if ((n[0] >= bias * n[1]) && (n[1] >= n[2]))
        {
            imax0 = 0;
            imax1 = 1;
            imax2 = 2;
        }

        // Try first with the largest dimension
        if (!div)
            for (int k = pmax - 1; k >= 0; k--)
            {
                if (fac[imax0][k] > 0 && fac[3][k] > 0)
                {
#ifndef NDEBUG
                    if (onpe0_ && os_ != 0)
                        (*os_) << "Divide domain along direction " << imax0
                               << endl;
#endif
                    n_mpi_tasks_dir_[imax0] *= coeff[k];
                    fac[imax0][k]--;
                    n[imax0] /= coeff[k];
                    fac[3][k]--;
                    n[3] /= coeff[k];
                    div = true;
                    break;
                }
            }
        if (!div)
            for (int k = pmax - 1; k >= 0; k--)
            {
                if (fac[imax1][k] > 0 && fac[3][k] > 0)
                {
#ifndef NDEBUG
                    if (onpe0_ && os_ != 0)
                        (*os_) << "Divide domain along direction " << imax1
                               << endl;
#endif
                    n_mpi_tasks_dir_[imax1] *= coeff[k];
                    fac[imax1][k]--;
                    n[imax1] /= coeff[k];
                    fac[3][k]--;
                    n[3] /= coeff[k];
                    div = true;
                    break;
                }
            }
        if (!div)
            for (int k = pmax - 1; k >= 0; k--)
            {
                if (fac[imax2][k] > 0 && fac[3][k] > 0)
                {
#ifndef NDEBUG
                    if (onpe0_ && os_ != 0)
                        (*os_) << "Divide domain along direction " << imax2
                               << endl;
#endif
                    n_mpi_tasks_dir_[imax2] *= coeff[k];
                    fac[imax2][k]--;
                    n[imax2] /= coeff[k];
                    fac[3][k]--;
                    n[3] /= coeff[k];
                    div = true;
                    break;
                }
            }
    }

#ifndef NDEBUG
    if (onpe0_ && os_ != 0)
    {
        (*os_) << " Subdomain: " << n[0] << " * " << n[1] << " * " << n[2]
               << endl;
        (*os_) << " MPI: " << n_mpi_tasks_dir_[0] << " * "
               << n_mpi_tasks_dir_[1] << " * " << n_mpi_tasks_dir_[2] << endl;
    }
#endif

    return n_mpi_tasks_dir_[0] * n_mpi_tasks_dir_[1] * n_mpi_tasks_dir_[2];
}

void PEenv::split_comm(const int nx, const int ny, const int nz, const int bias)
{
    int nmpi     = 1;
    int flag     = 1;
    color_       = 0;
    comm_active_ = comm_;

    MPI_Comm_rank(comm_, &mytask_);
    MPI_Comm_size(comm_, &n_mpi_tasks_);

    while (flag)
    {
        if (comm_active_ != comm_)
        {
            assert(comm_active_ != NULL);
            MPI_Comm_free(&comm_active_);
        }
        for (int i = 0; i < 3; i++)
        {
            n_mpi_tasks_dir_[i] = 1;
        }
        if (color_ == 0)
        {
            nmpi = geom(nx, ny, nz, bias);
            if (nmpi < n_mpi_tasks_
                && nmpi > 0) // reduces communicator size by 1
            {
                (*os_) << "WARNING!!! reduces communicator size by 1" << endl;
                if (mytask_ == (n_mpi_tasks_ - 1)) color_ = 1;
            }
        }
        int mpirc = MPI_Comm_split(comm_, color_, mytask_, &comm_active_);
        //(*os_)<<"MPI_Comm_split: comm_="<<comm_<<", color_="<<color_
        //    <<", mytask_="<<mytask_
        //    <<", comm_active_="<<comm_active_<<endl;
        if (mpirc != MPI_SUCCESS)
        {
            cerr << "MPI_Comm_split failed!, my color_=" << color_ << endl;
            MPI_Abort(comm_, 0);
        }
        MPI_Comm_size(comm_active_, &n_mpi_tasks_);
#ifndef NDEBUG
        if (color_ == 0 && onpe0_ && os_ != 0)
        {
            (*os_) << n_mpi_tasks_ << " MPI tasks" << endl;
        }
#endif

        // set flag to 0 once the number of processors matches nx*ny*nz
        if (nmpi == n_mpi_tasks_) flag = 0;
        if (nmpi == 0)
        {
            flag         = 0;
            n_mpi_tasks_ = 0;
            color_       = 1;
        }

        MPI_Bcast(&flag, 1, MPI_INT, 0, comm_);
    }
}

// Print list of node names
void PEenv::printPEnames(ostream& os) const
{
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    char buf[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    PMPI_Get_processor_name(processor_name, &namelen);
    if (mytask_ == 0)
        os << " Process " << mytask_ << " on " << processor_name << endl;

    if (color_ == 0)
        for (int ip = 1; ip < n_mpi_tasks_; ip++)
        {
            MPI_Barrier(comm_active_);
            if (mytask_ == 0)
            {
                MPI_Status status;
                int mpirc = MPI_Recv(&buf[0], MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
                    ip, ip, comm_active_, &status);
                if (mpirc != MPI_SUCCESS)
                {
                    cerr << "PEenv::printPEnames, MPI_Recv() failed!!!" << endl;
                    MPI_Abort(comm_, 0);
                }
            }
            else if (ip == mytask_)
            {
                // send processor name to pe0
                int mpirc = MPI_Send(&processor_name[0], MPI_MAX_PROCESSOR_NAME,
                    MPI_CHAR, 0, mytask_, comm_active_);
                if (mpirc != MPI_SUCCESS)
                {
                    cerr << "PEenv::printPEnames, MPI_Send() failed!!!" << endl;
                    MPI_Abort(comm_, 0);
                }
            }
            if (mytask_ == 0) os << " Process " << ip << " on " << buf << endl;
        }
    if (mytask_ == 0) os << endl;
}

void PEenv::globalExit(const int i) const
{
#ifdef USE_MPI
    MPI_Abort(comm_, 2);
#else
    exit(i);
#endif
}
#endif

void PEenv::bcast(int* val, const int n) const
{
#ifdef USE_MPI
    MPI_Bcast(val, n, MPI_INT, 0, comm_active_);
#endif
}

void PEenv::bcast(char* val, const int n) const
{
#ifdef USE_MPI
    MPI_Bcast(val, n, MPI_CHAR, 0, comm_active_);
#endif
}

void PEenv::bcast(double* val, const int n) const
{
#ifdef USE_MPI
    MPI_Bcast(val, n, MPI_DOUBLE, 0, comm_active_);
#endif
}

} // namespace pb
