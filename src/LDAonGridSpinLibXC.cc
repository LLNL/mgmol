// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifdef USE_LIBXC

#include "LDAonGridSpinLibXC.h"

#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Potentials.h"
#include "mputils.h"

template <class T>
void LDAonGridSpinLibXC<T>::update()
{
    get_xc_tm_.start();

    int iterative_index = rho_.getIterativeIndex();

    std::vector<double> rho(2 * np_);
    for (int j = 0; j < np_; j++)
    {
        for (short is = 0; is < 2; is++)
        {
            rho[2 * j + is] = rho_.rho_[is][j];
        }
    }

    std::vector<double> vxc(2 * np_);
    xc_lda_exc_vxc(&xfunc_, np_, rho.data(), &exc_[0], &vxc[0]);
    std::vector<double> vtmp(2 * np_);
    std::vector<double> etmp(np_);
    xc_lda_exc_vxc(&cfunc_, np_, rho.data(), &etmp[0], &vtmp[0]);

    int ione   = 1;
    double one = 1.;
    DAXPY(&np_, &one, &etmp[0], &ione, &exc_[0], &ione);

    int nn = np_ * 2;
    DAXPY(&nn, &one, &vxc[0], &ione, &vtmp[0], &ione);
    for (int j = 0; j < np_; j++)
    {
        vxc[j] = vtmp[2 * j + myspin_];
    }

    pot_.setVxc(&vxc[0], iterative_index);

    get_xc_tm_.stop();
}

template <class T>
double LDAonGridSpinLibXC<T>::getExc() const // in [Ha]
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    double sum = 0.;
    double exc = LinearAlgebraUtils<MemorySpace::Host>::MPdot(
        np_, &rho_.rho_[0][0], &exc_[0]);
    exc += LinearAlgebraUtils<MemorySpace::Host>::MPdot(
        np_, &rho_.rho_[1][0], &exc_[0]);
    exc *= mygrid.vel();

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int rc          = mmpi.allreduce(&exc, &sum, 1, MPI_SUM);
    if (rc != MPI_SUCCESS)
    {
        (*MPIdata::sout) << "MPI_Allreduce double sum failed!!!" << std::endl;
        mmpi.abort();
    }
    return sum;
}

template class LDAonGridSpinLibXC<LocGridOrbitals>;
template class LDAonGridSpinLibXC<ExtendedGridOrbitals>;

#endif
