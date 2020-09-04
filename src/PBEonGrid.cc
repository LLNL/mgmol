// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "PBEonGrid.h"

#include "Control.h"
#include "Delh4.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "Mesh.h"
#include "PBh4.h"

#include "Potentials.h"

Timer XConGrid::get_xc_tm_("XConGrid::get_xc");

template <class T>
void PBEonGrid<T>::update()
{
    get_xc_tm_.start();

    std::vector<std::vector<RHODTYPE>>& vrho = rho_.rho_;
    //    int     ione=1;
    double one = 1.;

    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    pb::Grid newGrid(mygrid, 2);

    pb::FDoper<RHODTYPE>* myoper_del[3];
    myoper_del[0] = new pb::Delxh4<RHODTYPE>(newGrid);
    myoper_del[1] = new pb::Delyh4<RHODTYPE>(newGrid);
    myoper_del[2] = new pb::Delzh4<RHODTYPE>(newGrid);

    pb::GridFunc<RHODTYPE> gf_rho(newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    gf_rho.assign(&vrho[0][0]);

    pb::GridFunc<RHODTYPE> gf_tmp(newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
#ifdef USE_LIBXC
    double* tmp   = new double[np_];
    double* sigma = new double[np_];
    memset(sigma, 0, np_ * sizeof(double));
#else
    RHODTYPE* tmp = new RHODTYPE[np_];
#endif

    // compute grad rho
    for (int i = 0; i < 3; i++)
    {
        myoper_del[i]->apply(gf_rho, gf_tmp);
        // convert gf_tmp back into double*
        gf_tmp.init_vect(tmp, 'd');
#ifdef USE_LIBXC
        for (int j = 0; j < np_; j++)
        {
            sigma[j] += (tmp[j] * tmp[j]);
        }
#else
        pbe_->setGradRho(i, tmp);
#endif
    }

    const int iterative_index = rho_.getIterativeIndex();

    delete[] tmp;

    for (int i = 0; i < 3; i++)
        delete myoper_del[i];

#ifdef USE_LIBXC
    xc_gga_exc_vxc(
        &xfunc_, np_, &rho_.rho_[0][0], sigma, &exc_[0], &vxc_[0], &vsigma_[0]);
    std::vector<double> vtmp(np_);
    std::vector<double> etmp(np_);
    std::vector<double> stmp(np_);
    xc_gga_exc_vxc(
        &cfunc_, np_, &rho_.rho_[0][0], sigma, &etmp[0], &vtmp[0], &stmp[0]);
    delete[] sigma;

    DAXPY(&np_, &one, &vtmp[0], &ione, &vxc_[0], &ione);
    DAXPY(&np_, &one, &etmp[0], &ione, &exc_[0], &ione);
    DAXPY(&np_, &one, &stmp[0], &ione, &vsigma_[0], &ione);

    vector<POTDTYPE> vstmp(np_);
    MPcpy(&vstmp[0], &vsigma_[0], np_);
    pb::GridFunc<POTDTYPE> gf_vsigma(
        newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    gf_vsigma.assign(&vstmp[0]);
    gf_vsigma *= 2.;
#else
    pbe_->computeXC();

    pb::GridFunc<POTDTYPE> gf_vsigma(
        newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    gf_vsigma.assign(pbe_->pvxc2_);
    gf_vsigma *= -1.;
#endif

    POTDTYPE* tmp_pot = new POTDTYPE[np_];
    pb::GridFunc<POTDTYPE> gf_vxc(newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    pb::DielFunc<POTDTYPE> diel_vxc2(gf_vsigma);
    pb::PBh4<POTDTYPE> myoper(newGrid, diel_vxc2);
    // convert gf_rho to POTDTYPE
    pb::GridFunc<POTDTYPE> gf_lhs(gf_rho);
    myoper.apply(gf_lhs, gf_vxc);
    // convert gf_vxc back into a double*
    gf_vxc.init_vect(tmp_pot, 'd');

#ifdef USE_LIBXC
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        &np_, &one, &vxc_[0], &ione, tmp_pot, &ione);
#else
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        np_, one, pbe_->pvxc1_, tmp_pot);
#endif

    pot_.setVxc(tmp_pot, iterative_index);

    delete[] tmp_pot;

    get_xc_tm_.stop();
}

template <class T>
double PBEonGrid<T>::getExc() const
{
    assert(pbe_ != nullptr);

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

#ifdef USE_LIBXC
    int ione = 1;
    double exc
        = mygrid.vel() * ddot(&np_, &rho_.rho_[0][0], &ione, &exc_[0], &ione);

    double sum      = 0.;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int rc          = mmpi.allreduce(&exc, &sum, 1, MPI_SUM);
    if (rc != MPI_SUCCESS)
    {
        (*MPIdata::sout) << "MPI_Allreduce double sum failed!!!" << endl;
        Control& ct = *(Control::instance());
        ct.global_exit(2);
    }
    return sum;
#else
    return mygrid.vel() * pbe_->computeRhoDotExc();
#endif
}

template class PBEonGrid<LocGridOrbitals>;
template class PBEonGrid<ExtendedGridOrbitals>;
