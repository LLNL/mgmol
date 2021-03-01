// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "PBEonGridSpin.h"

#include "Control.h"
#include "Delh4.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol_MPI.h"
#include "Mesh.h"
#include "PBh4.h"
#include "Potentials.h"
#include "mputils.h"

template <class T>
PBEonGridSpin<T>::PBEonGridSpin(Rho<T>& rho, Potentials& pot)
    : np_(rho.rho_[0].size()), rho_(rho), pot_(pot)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    myspin_         = mmpi.myspin();
    pbe_            = new PBEFunctional(rho.rho_);
    vxc_.resize(np_ * 2);
}

template <class T>
void PBEonGridSpin<T>::update()
{
    get_xc_tm_.start();

    int iterative_index = rho_.getIterativeIndex();

    std::vector<std::vector<RHODTYPE>>& vrho = rho_.rho_;

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

    RHODTYPE* grad_rho[2];
    for (short is = 0; is < 2; is++)
    {
        grad_rho[is] = new RHODTYPE[np_];
        memset(grad_rho[is], 0, np_ * sizeof(RHODTYPE));
    }

    pb::GridFunc<RHODTYPE> gf_tmp(newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    // compute grad rho, one direction at a time
    for (short dir = 0; dir < 3; dir++)
    {
        for (short is = 0; is < 2; is++)
        {
            memset(grad_rho[is], 0, np_ * sizeof(RHODTYPE));
            // gf_rho[is]->trade_boundaries();
            myoper_del[dir]->apply(*gf_rho[is], gf_tmp);
            // convert gf_tmp back into RHODTYPE*
            gf_tmp.init_vect(grad_rho[is], 'd');
        }
        pbe_->setGradRhoUp(dir, grad_rho[0]);
        pbe_->setGradRhoDn(dir, grad_rho[1]);
    }

    for (short is = 0; is < 2; is++)
        delete[] grad_rho[is];

    for (short i = 0; i < 3; i++)
        delete myoper_del[i];

    pb::GridFunc<POTDTYPE>* gf_vsigma[2];
    pbe_->computeXC();

    if (myspin_ == 0)
    {
        assert(pbe_->pvxc1_up_ != nullptr);
        MPcpy(&vxc_[0], pbe_->pvxc1_up_, np_);
        gf_vsigma[0] = new pb::GridFunc<POTDTYPE>(
            newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
        gf_vsigma[0]->assign(pbe_->pvxc2_upup_, 'd');

        gf_vsigma[1] = new pb::GridFunc<POTDTYPE>(
            newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
        gf_vsigma[1]->assign(pbe_->pvxc2_updn_, 'd');
    }
    else
    {
        assert(pbe_->pvxc1_dn_ != nullptr);
        MPcpy(&vxc_[np_], pbe_->pvxc1_dn_, np_);
        gf_vsigma[0] = new pb::GridFunc<POTDTYPE>(
            newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
        gf_vsigma[0]->assign(pbe_->pvxc2_dnup_, 'd');

        gf_vsigma[1] = new pb::GridFunc<POTDTYPE>(
            newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
        gf_vsigma[1]->assign(pbe_->pvxc2_dndn_, 'd');
    }
    (*gf_vsigma[0]) *= -1.;
    (*gf_vsigma[1]) *= -1.;

    POTDTYPE* tmp = new POTDTYPE[np_];
    pb::GridFunc<POTDTYPE> gf_tmp_pot(gf_tmp);
    for (short is = 0; is < 2; is++)
    {
        pb::DielFunc<POTDTYPE> diel(*gf_vsigma[is]);
        pb::PBh4<POTDTYPE> myoper(newGrid, diel);
        // convert gf_rho to POTDTYPE
        pb::GridFunc<POTDTYPE> gf_lhs(*gf_rho[is]);
        myoper.apply(gf_lhs, gf_tmp_pot);
        // convert gf_vxc back into a POTDTYPE*
        gf_tmp_pot.init_vect(tmp, 'd');
        LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
            np_, 1., tmp, &vxc_[np_ * myspin_]);
    }

    pot_.setVxc(&vxc_[np_ * myspin_], iterative_index);

    delete[] tmp;

    for (short isp = 0; isp < 2; isp++)
    {
        delete gf_rho[isp];
        delete gf_vsigma[isp];
    }
    get_xc_tm_.stop();
}

template <class T>
double PBEonGridSpin<T>::getExc() const
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    assert(pbe_ != nullptr);
    double exc = pbe_->computeRhoDotExc();
    return exc * mygrid.vel();
}

template class PBEonGridSpin<LocGridOrbitals>;
template class PBEonGridSpin<ExtendedGridOrbitals>;
