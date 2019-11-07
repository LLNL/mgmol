// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Preconditioning.h"
#include "LapFactory.h"

using namespace std;

template <typename T>
Preconditioning<T>::Preconditioning(const short lap_type, const short maxlevels,
    const pb::Grid& grid, const short bc[3])
{
    max_levels_ = maxlevels;
    lap_type_   = lap_type;
    for (short i = 0; i < 3; i++)
        bc_[i] = bc[i];

    pb::Grid* mygrid = new pb::Grid(grid);
    grid_.push_back(mygrid);

    pb::Lap<T>* myoper = LapFactory<T>::createLap(*grid_[0], lap_type);
    precond_oper_.push_back(myoper);
}

template <typename T>
Preconditioning<T>::Preconditioning(const Preconditioning& precond)
{
    max_levels_      = precond.max_levels_;
    lap_type_        = precond.lap_type_;
    pb::Grid* mygrid = new pb::Grid(*(precond.grid_[0]));
    grid_.push_back(mygrid);
    for (short i = 0; i < 3; i++)
        bc_[i] = precond.bc_[i];
}

template <typename T>
Preconditioning<T>::~Preconditioning()
{
    clear();
}

template <typename T>
void Preconditioning<T>::clear()
{
    for (short i = 0; i < (short)precond_oper_.size(); i++)
    {
        delete precond_oper_[i];
    }
    for (short i = 0; i < (short)gf_work_.size(); i++)
    {
        assert(gf_work_[i] != nullptr);
        delete gf_work_[i];
    }
    for (short i = 0; i < (short)gf_rcoarse_.size(); i++)
    {
        assert(gf_rcoarse_[i] != nullptr);
        delete gf_rcoarse_[i];
    }
    for (short i = 0; i < (short)gf_newv_.size(); i++)
    {
        assert(gf_newv_[i] != nullptr);
        delete gf_newv_[i];
    }
    // delete grids after pb::GridFunc<T> objects since those
    // have data members references to grids
    for (short i = 0; i < (short)grid_.size(); i++)
    {
        delete grid_[i];
    }
    precond_oper_.clear();
    grid_.clear();
    gf_work_.clear();
    gf_rcoarse_.clear();
    gf_newv_.clear();

    for (unsigned short l = 0; l < gfv_work_.size(); l++)
    {
        delete gfv_work_[l];
    }
    for (unsigned short l = 0; l < gfv_rcoarse_.size(); l++)
    {
        delete gfv_rcoarse_[l];
        delete gfv_newv_[l];
    }
    gfv_work_.clear();
    gfv_rcoarse_.clear();
    gfv_newv_.clear();
}

template <typename T>
void Preconditioning<T>::reset(map<int, GridMask*>& st_to_mask,
    const vector<vector<int>>& overlapping_gids)
{
    clear();

    setup(st_to_mask, overlapping_gids);
}

template <typename T>
void Preconditioning<T>::setup(map<int, GridMask*>& st_to_mask,
    const vector<vector<int>>& overlapping_gids)
{
    assert(overlapping_gids.size() > 0);

    overlapping_gids_ = overlapping_gids;

    // fine level
    pb::GridFunc<T>* gf_work
        = new pb::GridFunc<T>(*grid_[0], bc_[0], bc_[1], bc_[2]);
    gf_work_.push_back(gf_work);

    pb::GridFuncVector<T>* gfv_work = new pb::GridFuncVector<T>(
        true, *grid_[0], bc_[0], bc_[1], bc_[2], overlapping_gids_);
    gfv_work_.push_back(gfv_work);

    // coarse levels
    pb::Grid* mygrid = grid_[0];
    for (short ln = 1; ln <= max_levels_; ln++)
    {
        pb::Grid* coarse_grid = new pb::Grid(mygrid->coarse_grid());
        grid_.push_back(coarse_grid);

        // use 2nd order FD operator
        pb::Lap<T>* myoper = LapFactory<T>::createLap(*coarse_grid, 1);
        precond_oper_.push_back(myoper);

        gf_work = new pb::GridFunc<T>(*coarse_grid, bc_[0], bc_[1], bc_[2]);
        gf_work_.push_back(gf_work);

        pb::GridFunc<T>* gf_rcoarse
            = new pb::GridFunc<T>(*coarse_grid, bc_[0], bc_[1], bc_[2]);
        gf_rcoarse_.push_back(gf_rcoarse);
        pb::GridFunc<T>* gf_newv
            = new pb::GridFunc<T>(*coarse_grid, bc_[0], bc_[1], bc_[2]);
        gf_newv_.push_back(gf_newv);

        mygrid = coarse_grid;

        gfv_work = new pb::GridFuncVector<T>(
            true, *coarse_grid, bc_[0], bc_[1], bc_[2], overlapping_gids_);
        gfv_work_.push_back(gfv_work);

        pb::GridFuncVector<T>* gfv_rcoarse = new pb::GridFuncVector<T>(
            true, *coarse_grid, bc_[0], bc_[1], bc_[2], overlapping_gids_);
        gfv_rcoarse_.push_back(gfv_rcoarse);

        pb::GridFuncVector<T>* gfv_newv = new pb::GridFuncVector<T>(
            true, *coarse_grid, bc_[0], bc_[1], bc_[2], overlapping_gids_);
        gfv_newv_.push_back(gfv_newv);
    }
    gid2mask_ = st_to_mask;
}

// MG V-cycle with mask corresponding to state istate
// (no mask if istate==-1)
template <typename T>
void Preconditioning<T>::mg(pb::GridFuncVector<T>& gfv_v,
    const pb::GridFuncVector<T>& gfv_f, const short level)
{
#ifdef PRINT_OPERATIONS
    if (onpe0)
        (*MPIdata::sout) << "Preconditioning::mg() at level " << level << endl;
#endif
    assert(gfv_work_.size() > level);
    assert(gfv_work_[level] != nullptr);

    short ncycl = 2;
    if (level == max_levels_) ncycl = 4;

    pb::Lap<T>* myoper = precond_oper_[level];

    // SMOOTHING
    for (short it = 0; it < ncycl; it++)
    {
        myoper->jacobi(gfv_v, gfv_f, *gfv_work_[level]);
        app_mask(gfv_v, level);
    }

    if (level == max_levels_) return;

    // COARSE GRID CORRECTION

    // LOCALIZATION
    app_mask(*gfv_work_[level], level);

    // restrictions
    pb::GridFuncVector<T>* rcoarse = gfv_rcoarse_[level];
    assert(rcoarse != nullptr);
    gfv_work_[level]->restrict3D(*rcoarse);

    // LOCALIZATION
    app_mask(*rcoarse, level + 1);

    // storage functions for coarse grid
    pb::GridFuncVector<T>* newv = gfv_newv_[level];

    // call mgrid solver on a coarser level
    newv->resetData();
    mg(*newv, *rcoarse, level + 1);

    gfv_work_[level]->extend3D(*newv);

    // LOCALIZATION
    app_mask(*gfv_work_[level], level);

    gfv_v -= (*gfv_work_[level]);

    // post-smoothing
    for (short it = 0; it < 2; it++)
    {
        myoper->jacobi(gfv_v, gfv_f, *gfv_work_[level]);
        app_mask(gfv_v, level);
    }

    if (bc_[0] != 1 || bc_[2] != 1 || bc_[2] != 1) gfv_v.trade_boundaries();
}

template <typename T>
void Preconditioning<T>::mg(pb::GridFunc<T>& gf_v, const pb::GridFunc<T>& gf_f,
    const short level, const int istate)
{
    //(*MPIdata::sout)<<"Preconditioning::mg() at level "<<level<<endl;
    short ncycl = 2;
    if (level == max_levels_) ncycl = 4;

    pb::Lap<T>* myoper = precond_oper_[level];

    // SMOOTHING
    for (short it = 0; it < ncycl; it++)
    {
        myoper->jacobi(gf_v, gf_f, *gf_work_[level]);
        app_mask(gf_v, level, istate);
    }

    if (level == max_levels_) return;

    // COARSE GRID CORRECTION

    // LOCALIZATION
    app_mask(*gf_work_[level], level, istate);

    // restrictions
    pb::GridFunc<T>* rcoarse = gf_rcoarse_[level];
    gf_work_[level]->restrict3D(*rcoarse);

    // LOCALIZATION
    app_mask(*rcoarse, level + 1, istate);

    // storage functions for coarse grid
    pb::GridFunc<T>* newv = gf_newv_[level];

    // call mgrid solver on a coarser level
    newv->resetData();
    mg(*newv, *rcoarse, level + 1, istate);

    gf_work_[level]->extend3D(*newv);

    // LOCALIZATION
    app_mask(*gf_work_[level], level, istate);

    gf_v -= (*gf_work_[level]);

    // post-smoothing
    for (short it = 0; it < 2; it++)
    {
        myoper->jacobi(gf_v, gf_f, *gf_work_[level]);
        app_mask(gf_v, level, istate);
    }

    if (bc_[0] != 1 || bc_[2] != 1 || bc_[2] != 1) gf_v.trade_boundaries();
}

template <typename T>
void Preconditioning<T>::app_mask(
    pb::GridFuncVector<T>& gvu, const short level) const
{
    const int nfunc = (int)gvu.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int k = 0; k < nfunc; k++)
    {
        app_mask(gvu.func(k), level, k);
    }
}

template <typename T>
void Preconditioning<T>::app_mask(
    pb::GridFunc<T>& gu, const short level, const int color) const
{
    if (color == -1) return; // no mask applied
    if (gid2mask_.empty()) return; // no mask applied

    assert(overlapping_gids_.size() > 0);

    const short subdivx = overlapping_gids_.size();

    const pb::Grid& lgrid(gu.grid());
    const short nghosts = lgrid.ghost_pt();
    const int dim0      = lgrid.dim(0) / subdivx;
    const int incx      = lgrid.inc(0);
    const int lnumpt    = dim0 * incx;

    for (short iloc = 0; iloc < subdivx; iloc++)
    {
        int st = overlapping_gids_[iloc][color];

        if (st != -1)
        {
            map<int, GridMask*>::const_iterator it = gid2mask_.find(st);
            assert(it != gid2mask_.end());
            (it->second)->apply(gu, level, iloc);
        }
        else
        {
            int offset = (nghosts + dim0 * iloc) * incx;
            assert(offset + lnumpt < static_cast<int>(lgrid.sizeg()));
            T* pu = gu.uu() + offset;
            memset(pu, 0, lnumpt * sizeof(T));
        }
    }
}

template class Preconditioning<float>;
// template class Preconditioning<double>;
