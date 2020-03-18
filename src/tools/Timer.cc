// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Timer.h"

double Timer::gtod(void) const
{
    struct timeval tv;
    gettimeofday(&tv, (struct timezone*)nullptr);
    return tv.tv_sec + 1.e-6 * tv.tv_usec;
}

void Timer::print(std::ostream& os) const
{
    double tmin, tmax, tavg;

    double treal = total_real_;

    // get rank
    int mype;
    MPI_Comm_rank(comm_, &mype);

    int npes;
    MPI_Comm_size(comm_, &npes);

    // get min and max real time
    MPI_Allreduce(&treal, &tmin, 1, MPI_DOUBLE, MPI_MIN, comm_);
    MPI_Allreduce(&treal, &tmax, 1, MPI_DOUBLE, MPI_MAX, comm_);
    MPI_Allreduce(&treal, &tavg, 1, MPI_DOUBLE, MPI_SUM, comm_);
    tavg /= (double)(npes);

    // get min and max number of calls
    int numcalls = ncalls_;
    int nmin, nmax, navg;
    MPI_Allreduce(&numcalls, &nmin, 1, MPI_INT, MPI_MIN, comm_);
    MPI_Allreduce(&numcalls, &nmax, 1, MPI_INT, MPI_MAX, comm_);
    MPI_Allreduce(&numcalls, &navg, 1, MPI_INT, MPI_SUM, comm_);
    navg = (double)(navg) / (double)(npes);

    if (nmax == 0) return;

    if (mype == 0)
    {
        os.setf(std::ios::left, std::ios::adjustfield);
        os << "Timer: " << std::setw(50) << name_ << std::scientific
           << std::setprecision(2) << std::setw(9) << tmin << " / "
           << std::setprecision(2) << std::setw(9) << tavg << " / "
           << std::setprecision(2) << std::setw(9) << tmax << " / "
           << std::setprecision(2) << std::setw(7) << nmin << " / "
           << std::setw(7) << navg << " / " << std::setw(7) << nmax
           << std::endl;
    }
}
