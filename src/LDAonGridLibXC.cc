// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifdef USE_LIBXC

#include "LDAonGridLibXC.h"

#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Potentials.h"

template <class T>
void LDAonGridLibXC<T>::update()
{
    get_xc_tm_.start();

    int iterative_index = rho_.getIterativeIndex();

    xc_lda_exc_vxc(&xfunc_, exc_.size(), &rho_.rho_[0][0], &exc_[0], &vxc_[0]);
    std::vector<double> vtmp(vxc_.size());
    std::vector<double> etmp(vxc_.size());
    xc_lda_exc_vxc(&cfunc_, exc_.size(), &rho_.rho_[0][0], &etmp[0], &vtmp[0]);

    int ione   = 1;
    double one = 1.;
    int np     = vxc_.size();
    DAXPY(&np, &one, &vtmp[0], &ione, &vxc_[0], &ione);
    DAXPY(&np, &one, &etmp[0], &ione, &exc_[0], &ione);

    pot_.setVxc(&vxc_[0], iterative_index);

    get_xc_tm_.stop();
}

template <class T>
double LDAonGridLibXC<T>::getExc() const
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    int np = exc_.size();
    //        int ione=1;
    double exc = mygrid.vel()
                 * LinearAlgebraUtils<MemorySpace::Host>::MPdot(
                     np, &rho_.rho_[0][0], &exc_[0]);

    double sum      = 0.;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int rc          = mmpi.allreduce(&exc, &sum, 1, MPI_SUM);
    if (rc != MPI_SUCCESS)
    {
        (*MPIdata::sout) << "MPI_Allreduce double sum failed!!!" << std::endl;
        mmpi.abort();
    }
    return sum;
}

template class LDAonGridLibXC<LocGridOrbitals>;
template class LDAonGridLibXC<ExtendedGridOrbitals>;

#endif
