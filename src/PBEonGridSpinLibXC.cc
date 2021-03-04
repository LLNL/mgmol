// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifdef USE_LIBXC

#include "PBEonGridSpinLibXC.h"

#include "Control.h"
#include "Delh4.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol_MPI.h"
#include "Mesh.h"
#include "PBh4.h"
#include "Potentials.h"
#include "mputils.h"

#define WHITEBIRD94 1

template <class T>
PBEonGridSpinLibXC<T>::PBEonGridSpinLibXC(Rho<T>& rho, Potentials& pot)
    : np_(rho.rho_[0].size()), rho_(rho), pot_(pot)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    myspin_         = mmpi.myspin();
    int func_id     = XC_GGA_X_PBE;
    if (xc_func_init(&xfunc_, func_id, XC_POLARIZED) != 0)
    {
        std::cerr << "Functional " << func_id << " not found" << std::endl;
    }
    func_id = XC_GGA_C_PBE;
    if (xc_func_init(&cfunc_, func_id, XC_POLARIZED) != 0)
    {
        std::cerr << "Functional " << func_id << " not found" << std::endl;
    }
    exc_.resize(np_);
    vsigma_.resize(np_ * 3);
    vxc_.resize(np_ * 2);
}

template <class T>
void PBEonGridSpinLibXC<T>::update()
{
    get_xc_tm_.start();

    int iterative_index = rho_.getIterativeIndex();

    std::vector<std::vector<RHODTYPE>>& vrho = rho_.rho_;
    double one                               = 1.;

    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    pb::Grid newGrid(mygrid, 2);

    pb::FDoper<RHODTYPE>* myoper_del[3];
    myoper_del[0] = new pb::Delxh4<RHODTYPE>(newGrid);
    myoper_del[1] = new pb::Delyh4<RHODTYPE>(newGrid);
    myoper_del[2] = new pb::Delzh4<RHODTYPE>(newGrid);

    pb::GridFunc<RHODTYPE>* gf_rho[2];
    for (short is = 0; is < 2; is++)
    {
        gf_rho[is] = new pb::GridFunc<RHODTYPE>(
            newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
        gf_rho[is]->assign(&vrho[is][0], 'd');
    }

    std::vector<std::vector<std::vector<RHODTYPE>>> grad_rho;
    grad_rho.resize(2);
    for (short is = 0; is < 2; is++)
    {
        grad_rho[is].resize(3);
        for (short dir = 0; dir < 3; dir++)
            grad_rho[is][dir].resize(np_);
    }

    std::vector<double> sigma(3 * np_);
    memset(sigma.data(), 0, 3 * np_ * sizeof(double));
    std::vector<double> rho(2 * np_);
    for (int j = 0; j < np_; j++)
    {
        for (short is = 0; is < 2; is++)
        {
            rho[2 * j + is] = rho_.rho_[is][j];
        }
    }

    pb::GridFunc<RHODTYPE> gf_tmp(newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    // compute grad rho, one direction at a time
    for (short dir = 0; dir < 3; dir++)
    {
        for (short is = 0; is < 2; is++)
        {
            // gf_rho[is]->trade_boundaries();
            myoper_del[dir]->apply(*gf_rho[is], gf_tmp);
            // convert gf_tmp back into RHODTYPE*
            gf_tmp.init_vect(grad_rho[is][dir].data(), 'd');
        }
    }
    for (short dir = 0; dir < 3; dir++)
    {
        for (short is1 = 0; is1 < 2; is1++)
        {
            for (short is2 = is1; is2 < 2; is2++)
            {
                int jj = is1 + is2;
                for (int j = 0; j < np_; j++)
                {
                    sigma[jj]
                        += (grad_rho[is1][dir][j] * grad_rho[is2][dir][j]);
                    jj += 3;
                }
            }
        }
    }

    pb::GridFunc<POTDTYPE>* gf_vsigma[2];
    std::vector<double> vtmp(2 * np_);
    std::vector<double> etmp(np_);
    std::vector<double> stmp(3 * np_);
    xc_gga_exc_vxc(&xfunc_, np_, rho.data(), sigma.data(), &exc_[0], &vtmp[0],
        &vsigma_[0]);
    xc_gga_exc_vxc(
        &cfunc_, np_, rho.data(), sigma.data(), &etmp[0], &vxc_[0], &stmp[0]);

    int nn   = np_ * 2;
    int ione = 1;
    DAXPY(&nn, &one, &vxc_[0], &ione, &vtmp[0], &ione);
    nn = np_;
    DAXPY(&nn, &one, &etmp[0], &ione, &exc_[0], &ione);

    nn = np_ * 3;
    DAXPY(&nn, &one, &stmp[0], &ione, &vsigma_[0], &ione);

    // factor 2.
    double two = 2.;
    int ithree = 3;
    DSCAL(&np_, &two, &vsigma_[0], &ithree);
    // DSCAL(&np_, &two, &vsigma_[np_], &ithree);
    DSCAL(&np_, &two, &vsigma_[2], &ithree);

    if (myspin_ == 0)
    {
        for (int j = 0; j < np_; j++)
        {
            vxc_[j] = vtmp[2 * j];
        }
    }
    else
    {
        for (int j = 0; j < np_; j++)
        {
            vxc_[np_ + j] = vtmp[2 * j + 1];
        }
    }
    std::vector<POTDTYPE> vstmp(np_);
    gf_vsigma[0] = new pb::GridFunc<POTDTYPE>(
        newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    gf_vsigma[1] = new pb::GridFunc<POTDTYPE>(
        newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);

    if (myspin_ == 0)
    {
        for (int j = 0; j < np_; j++)
        {
            vstmp[j] = vsigma_[3 * j];
        }
        gf_vsigma[0]->assign(&vstmp[0], 'd');
        for (int j = 0; j < np_; j++)
        {
            vstmp[j] = vsigma_[3 * j + 1];
        }
        gf_vsigma[1]->assign(&vstmp[0], 'd');
    }
    else
    {
        for (int j = 0; j < np_; j++)
        {
            vstmp[j] = vsigma_[3 * j + 1];
        }
        gf_vsigma[0]->assign(&vstmp[0], 'd');

        for (int j = 0; j < np_; j++)
        {
            vstmp[j] = vsigma_[3 * j + 2];
        }
        gf_vsigma[1]->assign(&vstmp[0], 'd');
    }

    std::vector<POTDTYPE> tmp(np_);
#ifdef WHITEBIRD94
    //
    // White & Bird (1994)
    //
    std::vector<double> vsgrad(np_);
    pb::GridFunc<POTDTYPE> gf_vsgrad(
        newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    for (short is = 0; is < 2; is++)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < np_; j++)
            {
                vsgrad[j] = vsigma_[3 * j + is + myspin_] * grad_rho[is][i][j];
            }
            gf_vsgrad.assign(vsgrad.data());
            myoper_del[i]->apply(gf_vsgrad, gf_tmp);

            gf_tmp.init_vect(tmp.data(), 'd');
            LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
                np_, -1., tmp.data(), &vxc_[np_ * myspin_]);
        }
    }
#else
    pb::GridFunc<POTDTYPE> gf_tmp_pot(gf_tmp);
    for (short is = 0; is < 2; is++)
    {
        pb::DielFunc<POTDTYPE> diel(*gf_vsigma[is]);
        pb::PBh4<POTDTYPE> myoper(newGrid, diel);
        // convert gf_rho to POTDTYPE
        pb::GridFunc<POTDTYPE> gf_lhs(*gf_rho[is]);
        myoper.apply(gf_lhs, gf_tmp_pot);
        // convert gf_vxc back into a POTDTYPE*
        gf_tmp_pot.init_vect(tmp.data(), 'd');
        LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
            np_, one, tmp.data(), &vxc_[np_ * myspin_]);
    }
#endif

    pot_.setVxc(&vxc_[np_ * myspin_], iterative_index);

    for (short i = 0; i < 3; i++)
        delete myoper_del[i];

    for (short isp = 0; isp < 2; isp++)
    {
        delete gf_rho[isp];
        delete gf_vsigma[isp];
    }

    get_xc_tm_.stop();
}

template <class T>
double PBEonGridSpinLibXC<T>::getExc() const
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    double sum = 0.;
    double exc = LinearAlgebraUtils<MemorySpace::Host>::MPdot(
        np_, &rho_.rho_[0][0], &exc_[0]);
    exc += LinearAlgebraUtils<MemorySpace::Host>::MPdot(
        np_, &rho_.rho_[1][0], &exc_[0]);
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&exc, &sum, 1, MPI_SUM);
    exc = sum;
    return exc * mygrid.vel();
}

template class PBEonGridSpinLibXC<LocGridOrbitals>;
template class PBEonGridSpinLibXC<ExtendedGridOrbitals>;

#endif
