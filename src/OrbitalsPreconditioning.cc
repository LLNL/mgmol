// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "OrbitalsPreconditioning.h"

#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "MasksSet.h"
#include "Mesh.h"
#include "Potentials.h"
#include "Preconditioning.h"
#include "ProjectedMatricesInterface.h"


template <class T>
OrbitalsPreconditioning<T>::~OrbitalsPreconditioning()
{
    assert(is_set_);
    assert(precond_ != 0);

    delete precond_;
    if (gfv_work_ != 0)
    {
        delete gfv_work_;
        gfv_work_ = 0;
    }
    if (data_wghosts_ != 0 && mixed_precision_)
    {
        delete data_wghosts_;
        data_wghosts_ = 0;
    }
}

template <class T>
void OrbitalsPreconditioning<T>::setup(T& orbitals,
    const short mg_levels, const short lap_type, MasksSet* currentMasks,
    LocalizationRegions* lrs)
{
    assert(!is_set_);

    Control& ct(*(Control::instance()));
    Mesh* mymesh = Mesh::instance();
    const pb::Grid& mygrid(mymesh->grid());

    precond_ = new Preconditioning<MGPRECONDTYPE>(
        lap_type, mg_levels, mygrid, ct.bc);

    std::map<int, GridMask*> gid_to_mask;
    const vector<int>& overlap_gids(lrs->getOverlapGids());
    for (vector<int>::const_iterator it = overlap_gids.begin();
         it != overlap_gids.end(); it++)
    {
        int gid         = (*it);
        GridMask* maski = currentMasks->get_pmask(gid);
        assert(maski != 0);
        gid_to_mask.insert(pair<int, GridMask*>(gid, maski));
    }
    precond_->setup(gid_to_mask, orbitals.getOverlappingGids());

    assert(
        orbitals.chromatic_number() == orbitals.getOverlappingGids()[0].size());

    if (ct.blockPrecond())
    {
        gfv_work_ = new pb::GridFuncVector<MGPRECONDTYPE>(true, mygrid,
            ct.bc[0], ct.bc[1], ct.bc[2], orbitals.getOverlappingGids());
    }

    if (mixed_precision_)
    {
        data_wghosts_ = new pb::GridFuncVector<MGPRECONDTYPE>(true, mygrid,
            ct.bc[0], ct.bc[1], ct.bc[2], orbitals.getOverlappingGids());
    }
    else
    {
        // cast to please compiler in mixed precision case (i.e. when this line
        // is not actually executed)
        data_wghosts_ = dynamic_cast<pb::GridFuncVector<MGPRECONDTYPE>*>(
            orbitals.getPtDataWGhosts());
    }
    is_set_ = true;

    assert(data_wghosts_);
}

template <class T>
void OrbitalsPreconditioning<T>::precond_mg(T& orbitals)
{
    assert(is_set_);
    assert(precond_ != 0);
    assert(gamma_ > 0.);

#ifdef PRINT_OPERATIONS
    if (onpe0) (*MPIdata::sout) << "T::precond_mg()..." << endl;
#endif
    precond_tm_.start();

    // store residual in GridFuncVector<T> container
    // used for ghost values (no ghost values needed)
    if (mixed_precision_)
        orbitals.setDataWithGhosts(data_wghosts_);
    else
        orbitals.setDataWithGhosts();
    // trade_boundaries();

    Control& ct(*(Control::instance()));
    if (!ct.blockPrecond())
    {
        Mesh* mymesh = Mesh::instance();
        const pb::Grid& mygrid(mymesh->grid());
        pb::GridFunc<MGPRECONDTYPE> gf_work(
            mygrid, ct.bc[0], ct.bc[1], ct.bc[2]);
        for (int i = 0; i < orbitals.chromatic_number(); i++)
        {
            gf_work.resetData();

            gf_work.axpy(gamma_, data_wghosts_->func(i));
            precond_->mg(gf_work, data_wghosts_->func(i), 0, i);

            // gf_work.init_vect(psi(i),'d');
            orbitals.setPsi(gf_work, i);
        }
    }
    else
    {
        // block-implemented preconditioner
        assert(gfv_work_ != 0);

        gfv_work_->resetData();

        gfv_work_->axpy(gamma_, *data_wghosts_);
        precond_->mg(*gfv_work_, *data_wghosts_, 0);

        orbitals.setPsi(*gfv_work_);
    }

#ifdef PRINT_OPERATIONS
    if (onpe0) (*MPIdata::sout)
        << "OrbitalsPreconditioning<T>::precond_mg() done" << endl;
#endif
    precond_tm_.stop();
}

template <class T>
void OrbitalsPreconditioning<T>::setGamma(const pb::Lap<ORBDTYPE>& lapOper,
    const Potentials& pot, const short mg_levels,
    ProjectedMatricesInterface* proj_matrices)
{
    assert(precond_ != 0);
    assert(is_set_);

    const double small_eig = proj_matrices->getLowestEigenvalue();
    double diag            = lapOper.invDiagEl();
    double vmax            = pot.max();

    // diag * 4^{N_level+1}
    // gamma = inverse of the largest eigenvalue for the low frequency error
    gamma_ = diag;
    for (short ln = 0; ln <= mg_levels; ln++)
    {
        gamma_ *= 4.;
    }
    gamma_ = 1.0 / (2.0 / gamma_ + fabs(vmax - small_eig));
#ifdef DEBUG
    Control& ct(*(Control::instance()));
    if (onpe0 && ct.verbose > 2)
        (*MPIdata::sout) << " time step for low frequencies corrections = "
                         << gamma_ << endl;
#endif
}

template <class T>
void OrbitalsPreconditioning<T>::printTimers(ostream& os)
{
    precond_tm_.print(os);
}

template class OrbitalsPreconditioning<LocGridOrbitals>;
template class OrbitalsPreconditioning<ExtendedGridOrbitals>;
