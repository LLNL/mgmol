// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "LBFGS.h"
#include "Control.h"
#include "Electrostatic.h"
#include "MGmol.h"
#include "MasksSet.h"
#include "Mesh.h"
#include "ProjectedMatrices.h"

template <class OrbitalsType>
LBFGS<OrbitalsType>::LBFGS(OrbitalsType** orbitals, Ions& ions,
    Rho<OrbitalsType>& rho, ConstraintSet& constraints,
    std::shared_ptr<LocalizationRegions> lrs, ClusterOrbitals* local_cluster,
    MasksSet& masks, MasksSet& corrmasks, Electrostatic& electrostat,
    const double dt, MGmol<OrbitalsType>& strategy)
    : IonicAlgorithm<OrbitalsType>(
        orbitals, ions, rho, constraints, lrs, masks, strategy),
      orbitals_(orbitals),
      ions_(ions),
      rho_(rho),
      lrs_(lrs),
      local_cluster_(local_cluster),
      masks_(masks),
      corrmasks_(corrmasks),
      electrostat_(electrostat),
      mgmol_strategy_(strategy)
{
    setup(dt);
}

template <class OrbitalsType>
void LBFGS<OrbitalsType>::setup(const double dt)
{
    if (lrs_)
        ref_lrs_ = std::shared_ptr<LocalizationRegions>(
            new LocalizationRegions(*lrs_));
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    Control& ct            = *(Control::instance());

    etot_i_[0] = 10000.;
    etot_i_[1] = 10000.;
    etot_i_[2] = 10000.;

    stepper_ = new LBFGS_IonicStepper(dt, IonicAlgorithm<OrbitalsType>::atmove_,
        IonicAlgorithm<OrbitalsType>::tau0_,
        IonicAlgorithm<OrbitalsType>::taup_,
        IonicAlgorithm<OrbitalsType>::fion_, IonicAlgorithm<OrbitalsType>::gid_,
        5, &etot_i_[0]);
    IonicAlgorithm<OrbitalsType>::registerStepper(stepper_);

    if (lrs_)
    {
        ref_masks_ = std::shared_ptr<MasksSet>(
            new MasksSet(lrs_, false, ct.getMGlevels()));
        ref_corrmasks_ = std::shared_ptr<MasksSet>(new MasksSet(lrs_, true, 0));
    }

    vh_init_ = new pb::GridFunc<POTDTYPE>(electrostat_.getVh());

    MasksSet* ref_masks     = ref_masks_ ? ref_masks_.get() : nullptr;
    MasksSet* ref_corrmasks = ref_corrmasks_ ? ref_corrmasks_.get() : nullptr;
    ref_orbitals_           = std::shared_ptr<OrbitalsType>(
        new OrbitalsType("LBFGS_ref", mygrid, mymesh->subdivx(), ct.numst,
                      ct.bcWF, (*orbitals_)->getProjMatrices(), ref_lrs_, ref_masks,
                      ref_corrmasks, local_cluster_));

    ref_orbitals_->assign(**orbitals_);
}

template <class OrbitalsType>
LBFGS<OrbitalsType>::~LBFGS()
{
    delete vh_init_;
    delete stepper_;
}

template <class OrbitalsType>
void LBFGS<OrbitalsType>::reset(const double dt)
{
    delete vh_init_;
    delete stepper_;

    setup(dt);
}

template <class OrbitalsType>
int LBFGS<OrbitalsType>::quenchElectrons(const int itmax, double& etot)
{
    etot_i_[0] = etot_i_[1];
    etot_i_[1] = etot_i_[2];
    int ret = mgmol_strategy_.quench(**orbitals_, ions_, itmax, 0, etot_i_[2]);

    etot = etot_i_[2];
    return ret;
}

template <class OrbitalsType>
void LBFGS<OrbitalsType>::updateRefs()
{
    if (!stepper_->check_last_step_accepted())
    {
        if (onpe0)
            (*MPIdata::sout)
                << "LBFGS: Reset orbitals to reference orbitals " << std::endl;
        (*orbitals_)->assign(*ref_orbitals_);
        electrostat_.setupInitialVh(*vh_init_);
    }
    else
    {
        // update references
        updateRefMasks();

        if (onpe0)
            (*MPIdata::sout)
                << "LBFGS: Update reference orbitals " << std::endl;
        ref_orbitals_->assign(**orbitals_);
        vh_init_->assign(electrostat_.getVh(), 'd');
    }
}

template <>
void LBFGS<LocGridOrbitals<ORBDTYPE>>::updateRefMasks()
{
    Control& ct = *(Control::instance());

    if (ct.lr_update)
    {
        assert(ref_lrs_);

        *ref_lrs_       = *lrs_;
        *ref_masks_     = masks_;
        *ref_corrmasks_ = corrmasks_;
    }

    MasksSet* ref_masks     = ref_masks_ ? ref_masks_.get() : nullptr;
    MasksSet* ref_corrmasks = ref_corrmasks_ ? ref_corrmasks_.get() : nullptr;
    ref_orbitals_->reset(ref_masks, ref_corrmasks, lrs_);
}

template <>
void LBFGS<ExtendedGridOrbitals<ORBDTYPE>>::updateRefMasks()
{
}

template <class OrbitalsType>
void LBFGS<OrbitalsType>::setQuenchTol() const
{
    Control& ct = *(Control::instance());
    ct.conv_tol = 0.1 * stepper_->etol();
    if (onpe0)
        (*MPIdata::sout) << std::setprecision(12) << std::fixed
                         << "LBFGS: Set SC convergence criterion to "
                         << ct.conv_tol << std::endl;
}

template <class OrbitalsType>
void LBFGS<OrbitalsType>::updatePotAndMasks()
{
    if (!stepper_->check_last_step_accepted())
        electrostat_.setupInitialVh(*vh_init_);

    IonicAlgorithm<OrbitalsType>::updatePotAndMasks();
}

template <class OrbitalsType>
bool LBFGS<OrbitalsType>::lbfgsLastStepNotAccepted() const
{
    return !stepper_->check_last_step_accepted();
}

template class LBFGS<LocGridOrbitals<ORBDTYPE>>;
template class LBFGS<ExtendedGridOrbitals<ORBDTYPE>>;
