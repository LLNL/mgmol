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

template <class T>
void PBEonGrid<T>::update()
{
    get_xc_tm_.start();

    std::vector<std::vector<RHODTYPE>>& vrho = rho_.rho_;

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
    std::vector<RHODTYPE> tmp(np_);

    // compute grad rho
    for (int i = 0; i < 3; i++)
    {
        myoper_del[i]->apply(gf_rho, gf_tmp);
        // convert gf_tmp back into double*
        gf_tmp.init_vect(tmp.data(), 'd');
        pbe_->setGradRho(i, tmp.data());
    }

    const int iterative_index = rho_.getIterativeIndex();

    for (int i = 0; i < 3; i++)
        delete myoper_del[i];

    // compute vxc1 and vxc2
    pbe_->computeXC();

    pb::GridFunc<POTDTYPE> gf_vxc2(newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    gf_vxc2.assign(pbe_->pvxc2_);
    gf_vxc2 *= -1.;

    //
    // compute vxc from vxc1 and vxc2
    //
    pb::GridFunc<POTDTYPE> gf_vxc(newGrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    pb::DielFunc<POTDTYPE> diel_vxc2(gf_vxc2);
    // convert gf_rho to POTDTYPE
    pb::GridFunc<POTDTYPE> gf_lhs(gf_rho);

    // compute div ( vxc2 * grad(n) )
    pb::PBh4<POTDTYPE> myoper(newGrid, diel_vxc2);
    myoper.apply(gf_lhs, gf_vxc);

    // convert gf_vxc back into data without ghosts
    std::vector<POTDTYPE> tmp_vxc(np_);
    gf_vxc.init_vect(tmp_vxc.data(), 'd');

    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        np_, 1., pbe_->pvxc1_, tmp_vxc.data());

    pot_.setVxc(tmp_vxc.data(), iterative_index);

    get_xc_tm_.stop();
}

template <class T>
double PBEonGrid<T>::getExc() const
{
    assert(pbe_ != nullptr);

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    return mygrid.vel() * pbe_->computeRhoDotExc();
}

template class PBEonGrid<LocGridOrbitals>;
template class PBEonGrid<ExtendedGridOrbitals>;
