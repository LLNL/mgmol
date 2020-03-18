// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "DirectionalReduce.h"

#include <cassert>
#include <math.h>

DirectionalReduce::DirectionalReduce(MPI_Comm cart_comm, const int max_steps[3])
    : cart_comm_(cart_comm)
{
#ifndef NDEBUG
    int status;
    int ret = MPI_Topo_test(cart_comm_, &status);
    assert(ret == MPI_SUCCESS);
    assert(status == MPI_CART);
#endif

    int periodic[3];
    int coords[3];
    MPI_Cart_get(cart_comm_, 3, nprocs_, periodic, coords);

    /* compute left and right steps in xyz directions */
    computeNumSteps(max_steps);
}

DirectionalReduce::DirectionalReduce(
    MPI_Comm cart_comm, const double spread_radius, const double domain[3])
    : cart_comm_(cart_comm)
{
    assert(spread_radius > 0.);

#ifndef NDEBUG
    int status;
    int ret = MPI_Topo_test(cart_comm_, &status);
    assert(ret == MPI_SUCCESS);
    assert(status == MPI_CART);
#endif

    int periodic[3];
    int coords[3];
    MPI_Cart_get(cart_comm_, 3, nprocs_, periodic, coords);

    /* compute left and right steps in xyz directions */
    computeNumSteps(spread_radius, domain);
}

void DirectionalReduce::computeDirMax(const short dir, int data[2])
{
    setupPersistentRequests(dir);

    sbuf_[0] = data[0];
    sbuf_[1] = data[1];

    /* Step in the left direction */
    short step = 0;
    while (step < lstep_[dir])
    {
        /* Send and receive data */
        MPI_Start(&request_[0]); // recv
        MPI_Start(&request_[1]); // send
        MPI_Waitall(2, request_, MPI_STATUSES_IGNORE);
        /* get max. sizes */
        sbuf_[0] = sbuf_[0] > rbuf_[0] ? sbuf_[0] : rbuf_[0];
        sbuf_[1] = sbuf_[1] > rbuf_[1] ? sbuf_[1] : rbuf_[1];
        step++;
    }

    /* Step in the right direction */
    step = 0;
    while (step < rstep_[dir])
    {
        /* Send and receive data */
        MPI_Start(&request_[2]); // recv
        MPI_Start(&request_[3]); // send
        MPI_Waitall(2, &request_[2], MPI_STATUSES_IGNORE);
        /* get max. sizes */
        sbuf_[0] = sbuf_[0] > rbuf_[0] ? sbuf_[0] : rbuf_[0];
        sbuf_[1] = sbuf_[1] > rbuf_[1] ? sbuf_[1] : rbuf_[1];
        step++;
    }

    data[0] = sbuf_[0];
    data[1] = sbuf_[1];

    deletePersistentRequests();
}

void DirectionalReduce::setupPersistentRequests(const short dir)
{
    // initialize persisted requests
    request_[0] = MPI_REQUEST_NULL;
    request_[1] = MPI_REQUEST_NULL;
    request_[2] = MPI_REQUEST_NULL;
    request_[3] = MPI_REQUEST_NULL;

    int source, dest;
    // Get source and destination ID to send and recv data - left direction
    short disp = -1;
    MPI_Cart_shift(cart_comm_, dir, disp, &source, &dest);
    /* setup request for left direction */
    MPI_Recv_init(&rbuf_[0], 2, MPI_INT, source, 0, cart_comm_, &request_[0]);
    MPI_Send_init(&sbuf_[0], 2, MPI_INT, dest, 0, cart_comm_, &request_[1]);

    // setup request for right direction
    // just reverse source and dest (no need for mpi_cart_shift())
    MPI_Recv_init(&rbuf_[0], 2, MPI_INT, dest, 1, cart_comm_, &request_[2]);
    MPI_Send_init(&sbuf_[0], 2, MPI_INT, source, 1, cart_comm_, &request_[3]);
}

void DirectionalReduce::deletePersistentRequests()
{
    MPI_Request_free(&request_[0]);
    MPI_Request_free(&request_[1]);
    MPI_Request_free(&request_[2]);
    MPI_Request_free(&request_[3]);
}

void DirectionalReduce::computeNumSteps(const int max_steps[3])
{
    assert(nprocs_[0] >= 1);
    assert(nprocs_[1] >= 1);
    assert(nprocs_[2] >= 1);

    /* compute left and right steps in xyz directions */

    /* x-direction */
    lstep_[0] = std::min(max_steps[0], (nprocs_[0] - 1));
    rstep_[0] = std::min(lstep_[0], nprocs_[0] - lstep_[0] - 1);

    /* y-direction */
    lstep_[1] = std::min(max_steps[1], (nprocs_[1] - 1));
    rstep_[1] = std::min(lstep_[1], nprocs_[1] - lstep_[1] - 1);

    /* z-direction */
    lstep_[2] = std::min(max_steps[2], (nprocs_[2] - 1));
    rstep_[2] = std::min(lstep_[2], nprocs_[2] - lstep_[2] - 1);
}

void DirectionalReduce::computeNumSteps(
    const double spread_radius, const double domain[3])
{
    assert(nprocs_[0] >= 1);
    assert(nprocs_[1] >= 1);
    assert(nprocs_[2] >= 1);

    assert(spread_radius > 0.);

    assert(domain[0] > 0.);
    assert(domain[1] > 0.);
    assert(domain[2] > 0.);

    /* First compute processor width info */
    double width[3];
    width[0] = domain[0] / nprocs_[0];
    width[1] = domain[1] / nprocs_[1];
    width[2] = domain[2] / nprocs_[2];

    /* compute left and right steps in xyz directions */
    // x-direction
    lstep_[0]
        = std::min((int)(ceil(spread_radius / width[0])), (nprocs_[0] - 1));
    rstep_[0] = std::min(lstep_[0], nprocs_[0] - lstep_[0] - 1);

    // y-direction
    lstep_[1]
        = std::min((int)(ceil(spread_radius / width[1])), (nprocs_[1] - 1));
    rstep_[1] = std::min(lstep_[1], nprocs_[1] - lstep_[1] - 1);

    // z-direction
    lstep_[2]
        = std::min((int)(ceil(spread_radius / width[2])), (nprocs_[2] - 1));
    rstep_[2] = std::min(lstep_[2], nprocs_[2] - lstep_[2] - 1);
}

void DirectionalReduce::printStats(std::ostream& os)
{
    os << "x-direction." << std::endl;
    os << "lstep = " << lstep_[0] << ", rstep = " << rstep_[0] << std::endl;
    os << "y-direction." << std::endl;
    os << "lstep = " << lstep_[1] << ", rstep = " << rstep_[1] << std::endl;
    os << "z-direction." << std::endl;
    os << "lstep = " << lstep_[2] << ", rstep = " << rstep_[2] << std::endl;
}
