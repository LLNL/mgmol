// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PB_PEENV_H
#define PB_PEENV_H

#include <cassert>
#include <iostream>
#include <unistd.h>

#include <mpi.h>

namespace pb
{

// x direction
#define EAST 0
#define WEST 1

// y direction
#define NORTH 2
#define SOUTH 3

// z direction
#define UP 4
#define DOWN 5

class PEenv
{
private:
    int mytask_;
    int n_mpi_tasks_;
    int n_mpi_tasks_dir_[3]; // nb tasks in x,y,z directions
    int mytask_dir_[3]; // task index in x,y,z directions
    int* other_tasks_dir_[3];
    int mpi_neighbors_[6];

    // base communicator from which all others are built
    const MPI_Comm comm_;

    MPI_Comm comm_active_;
    MPI_Comm cart_comm_;

    // communicator to exchange data along x, y, or z direction
    MPI_Comm comm_x_;
    MPI_Comm comm_y_;
    MPI_Comm comm_z_;
    int root_x_;

    int color_;
    bool onpe0_;
    std::ostream* os_;

    void setup_my_neighbors();
    int geom(const int, const int, const int, const int);
    void set_other_tasks_dir();
    void split_comm(const int, const int, const int, const int);

public:
    bool onpe0() const { return onpe0_; }
    int mytask() const { return mytask_; }
    int n_mpi_tasks() const { return n_mpi_tasks_; }
    int n_mpi_task(const int i) const { return n_mpi_tasks_dir_[i]; }
    int my_mpi(const int i) const { return mytask_dir_[i]; }
    int mpi_neighbors(int i) const { return mpi_neighbors_[i]; }
    MPI_Comm comm() const { return comm_active_; }
    MPI_Comm cart_comm() const { return cart_comm_; }
    MPI_Comm comm_active() const { return comm_active_; }
    MPI_Comm comm_x() const { return comm_x_; }
    MPI_Comm comm_y() const { return comm_y_; }
    MPI_Comm comm_z() const { return comm_z_; }
    MPI_Comm comm_global() const { return comm_; }
    int root_x() const { return root_x_; };
    int color() const { return color_; }
    int other_tasks_dir(const int dir, const int itask) const
    {
        return other_tasks_dir_[dir][itask];
    }

    PEenv(MPI_Comm comm, std::ostream* os = nullptr); // default constructor
    explicit PEenv(MPI_Comm comm, const int nx, const int ny, const int nz,
        const int bias = 4, std::ostream* os = nullptr);
    // explicit PEenv(MPI_Comm comm, const int mytask,
    //               const int taskx, const int tasky, const int taskz, ostream*
    //               os=0);

    ~PEenv();

    void task2xyz();

    void barrier(void) const;
    void bcast(char*, const int n) const;
    void bcast(int*, const int n) const;
    void bcast(double*, const int n) const;

    double double_sum_all(double) const;
    double double_max_all(double) const;
    double double_min_all(double) const;
    int int_max_all(int);

    int xyz2task(const int x, const int y, const int z) const
    {
#if 1
        assert(cart_comm_ != 0);
        int rank;
        int coords[3] = { x, y, z };
        MPI_Cart_rank(cart_comm_, coords, &rank);
        return rank;
#else
        return ((x + n_mpi_task(0)) % n_mpi_task(0)) * n_mpi_task(1)
                   * n_mpi_task(2)
               + ((y + n_mpi_task(1)) % n_mpi_task(1)) * n_mpi_task(2)
               + ((z + n_mpi_task(2)) % n_mpi_task(2));
#endif
    }

    void printPEnames(std::ostream&) const;

    void globalExit() const;

    void Isend(double* buf, int sizeb, const short dst, MPI_Request* req) const
    {
        MPI_Isend(
            buf, sizeb, MPI_DOUBLE, mpi_neighbors(dst), 0, cart_comm_, req);
    }
    void Isend(float* buf, int sizeb, const short dst, MPI_Request* req) const
    {
        MPI_Isend(
            buf, sizeb, MPI_FLOAT, mpi_neighbors(dst), 0, cart_comm_, req);
    }
    void Isend(int* buf, int sizeb, const short dst, MPI_Request* req) const
    {
        MPI_Isend(buf, sizeb, MPI_INT, mpi_neighbors(dst), 0, cart_comm_, req);
    }

    void Irecv(double* buf, int sizeb, const short src, MPI_Request* req) const
    {
        MPI_Irecv(
            buf, sizeb, MPI_DOUBLE, mpi_neighbors(src), 0, cart_comm_, req);
    }
    void Irecv(float* buf, int sizeb, const short src, MPI_Request* req) const
    {
        MPI_Irecv(
            buf, sizeb, MPI_FLOAT, mpi_neighbors(src), 0, cart_comm_, req);
    }
    void Irecv(int* buf, int sizeb, const short src, MPI_Request* req) const
    {
        MPI_Irecv(buf, sizeb, MPI_INT, mpi_neighbors(src), 0, cart_comm_, req);
    }

    // functions to reduce an int in one direction
    void maxXdir(int* values, const int count) const
    {
        int* sendbuf = new int[count];
        for (int i = 0; i < count; i++)
            sendbuf[i] = values[i];
        int mpierr = MPI_Allreduce(
            sendbuf, &values[0], count, MPI_INT, MPI_MAX, comm_x_);
        if (mpierr != MPI_SUCCESS)
        {
            std::cerr << "ERROR in PEenv::maxXdir()" << std::endl;
            sleep(5);
            MPI_Abort(comm_, EXIT_FAILURE);
        }
        delete[] sendbuf;
    }
    void maxYdir(int* values, const int count) const
    {
        int* sendbuf = new int[count];
        for (int i = 0; i < count; i++)
            sendbuf[i] = values[i];
        int mpierr = MPI_Allreduce(
            sendbuf, &values[0], count, MPI_INT, MPI_MAX, comm_y_);
        if (mpierr != MPI_SUCCESS)
        {
            std::cerr << "ERROR in PEenv::maxYdir()" << std::endl;
            sleep(5);
            MPI_Abort(comm_, EXIT_FAILURE);
        }
        delete[] sendbuf;
    }
    void maxZdir(int* values, const int count) const
    {
        int* sendbuf = new int[count];
        for (int i = 0; i < count; i++)
            sendbuf[i] = values[i];
        int mpierr = MPI_Allreduce(
            sendbuf, &values[0], count, MPI_INT, MPI_MAX, comm_z_);
        if (mpierr != MPI_SUCCESS)
        {
            std::cerr << "ERROR in PEenv::maxZdir()" << std::endl;
            sleep(5);
            MPI_Abort(comm_, EXIT_FAILURE);
        }
        delete[] sendbuf;
    }
};

} // namespace pb

#endif
